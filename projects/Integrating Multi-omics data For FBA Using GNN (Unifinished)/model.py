import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import k_hop_subgraph, shortest_path



class CustomGNNLayer(MessagePassing):
    def __init__(self, in_dim, out_dim):
        super().__init__(aggr='add')
        self.psi = nn.Linear(2*in_dim, out_dim)
        self.phi = nn.Linear(in_dim + out_dim, out_dim)
        
    def forward(self, x, edge_index):
        return self.propagate(edge_index, x=x)
    
    def message(self, x_i, x_j):
        return F.relu(self.psi(torch.cat([x_i, x_j], dim=-1)))
    
    def update(self, aggr, x):
        return F.relu(self.phi(torch.cat([x, aggr], dim=-1)))

class FinalGNNModel(nn.Module):
    def __init__(self, in_dim, hidden_dim, num_layers, num_heads, max_distance):
        super().__init__()
        self.gnn_layers = nn.ModuleList([
            CustomGNNLayer(hidden_dim if i>0 else in_dim, hidden_dim)
            for i in range(num_layers)
        ])
        self.num_heads = num_heads
        self.max_distance = max_distance
        self.attention_weights = nn.ParameterList([
            nn.Parameter(torch.randn(2*hidden_dim + max_distance))
            for _ in range(num_heads)
        ])
        self.final_layer = nn.Linear(num_heads*hidden_dim, 1)
        
    def forward(self, x, edge_index, edge_type, precomputed_dist):
        # 1. GNN Layers
        h = [x]
        for layer in self.gnn_layers:
            h.append(layer(h[-1], edge_index))
        
        # 2. Generate node pairs
        all_pairs = []
        for u in range(x.size(0)):
            for v in range(x.size(0)):
                if u != v and precomputed_dist[u][v] <= self.max_distance:
                    for l in range(len(h)-1):
                        for l_prime in range(len(h)-1):
                            all_pairs.append((
                                h[l][u], 
                                h[l_prime][v], 
                                precomputed_dist[u][v]
                            ))
        
        # 3. Multi-head attention
        outputs = []
        for head in range(self.num_heads):
            attn_scores = []
            projected = []
            for (hu, hv, d) in all_pairs:
                # One-hot encode distance
                d_onehot = F.one_hot(torch.tensor(d-1), self.max_distance)
                feature = torch.cat([hu, hv, d_onehot.float()])
                
                # Compute attention score
                score = torch.dot(self.attention_weights[head], feature)
                attn_scores.append(score)
                projected.append(self.attention_weights[head] * feature)
            
            # Softmax and aggregate
            attn_weights = F.softmax(torch.stack(attn_scores), dim=0)
            head_output = sum(w*p for w,p in zip(attn_weights, projected))
            outputs.append(head_output)
        
        # 4. Final classification
        combined = torch.cat(outputs, dim=-1)
        return torch.sigmoid(self.final_layer(combined))

# Example usage
model = FinalGNNModel(
    in_dim=64,
    hidden_dim=64,
    num_layers=3,
    num_heads=2,
    max_distance=5
)

# Precompute shortest paths (simplified example)
precomputed_dist = torch.randint(1, 6, (100, 100))  # Mock distance matrix

# Mock data
x = torch.randn(100, 64)  # Node features
edge_index = torch.randint(0, 100, (2, 500))       # Edge connections
edge_type = torch.randint(0, 2, (500,))            # 0=conjunctive, 1=disjunctive

output = model(x, edge_index, edge_type, precomputed_dist)
print(output.shape)  # [100, 1]
