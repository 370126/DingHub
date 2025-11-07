import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, ax = plt.subplots(figsize=(16, 12))

# Function to draw layer box
def draw_layer(x, y, w, h, label, detail, color="lightblue"):
    ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", edgecolor="black", facecolor=color))
    ax.text(x + w/2, y + h*0.65, label, ha='center', va='center', fontsize=11, weight='bold')
    ax.text(x + w/2, y + h*0.35, detail, ha='center', va='center', fontsize=9, wrap=True)

# Function to draw arrow between layers
def connect_layers(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', lw=1.5))

# Draw each module (in English)
draw_layer(0, 9, 3, 1.8, "Transcriptomic Input", "G = # common genes\nFeature: MinMax-avg across samples\nShape: [G, 1]", "lightgreen")
draw_layer(0, 6, 3, 1.8, "Epigenetic Input", "G = # common genes\nFeature: from epi_lung.csv\nShape: [G, 1]", "lightgreen")
draw_layer(4, 7.5, 3, 2.2, "Dual-Modality Graph", "2G nodes: G transcriptomic + G epigenetic\nEdges from label='conjunctive'\nedge_index Shape: [2, E]", "lightyellow")
draw_layer(8, 8.5, 3, 1.5, "GCNConv Layer 1", "Input: [2G, 1] → Output: [2G, 32]\nGraph convolution + ReLU", "lightblue")
draw_layer(11.5, 8.5, 3, 1.5, "GCNConv Layer 2", "Input: [2G, 32] → Output: [2G, 32]\nFurther propagation", "lightblue")
draw_layer(15, 8.5, 3, 1.5, "Attention Heads ×2", "Multi-head attention using edge distance\nOutput: [2G, 64]", "lightskyblue")
draw_layer(18.5, 8.5, 3, 1.5, "Modality Fusion + MLP", "Average transcriptomic & epigenetic node\nProject to expression ∈ [0, 10]\nShape: [G, 1]", "orange")
draw_layer(22, 8.5, 3, 1.5, "Differentiable GPR Bounds", "Parse COBRA GPR rules + expr.\nCompute lb / ub\nShape: [#Reactions]", "orchid")
draw_layer(25.5, 8.5, 3, 1.5, "FBA Layer (cvxpylayers)", "Input: bounds lb, ub\nOutput: predicted flux\nShape: [#Reactions]", "salmon")




# Draw arrows
connect_layers(3, 10, 4, 9.3)  # Transcriptomic → Graph
connect_layers(3, 7, 4, 8.2)   # Epigenetic → Graph
connect_layers(7, 8.6, 8, 9.1) # Graph → GCN1
connect_layers(11, 9.1, 11.5, 9.1) # GCN1 → GCN2
connect_layers(14.5, 9.1, 15, 9.1) # GCN2 → Attention
connect_layers(18, 9.1, 18.5, 9.1) # Attention → Projection
connect_layers(21.5, 9.1, 22, 9.1) # Projection → GPR
connect_layers(25, 9.1, 25.5, 9.1) # GPR → FBA

# Final flux output label
ax.text(29, 9.1, "Predicted\nFluxes", ha="center", va="center", fontsize=10, weight="bold", bbox=dict(boxstyle="round", fc="lightcoral", ec="black"))

# Style
ax.set_xlim(-1, 31)
ax.set_ylim(5, 11.2)
ax.axis('off')
plt.title("Multi-Omics GNN Framework for Flux Prediction", fontsize=15, weight='bold')
plt.tight_layout()
plt.show()

# save
fig.savefig('multi_omics_gnn_framework.png', dpi=300, bbox_inches='tight')




import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

# Function to draw boxes
def draw_box(ax, x, y, width, height, text, color='lightblue', fill=True):
    box = Rectangle((x - width/2, y - height/2), width, height, fill=fill, 
                    color=color, edgecolor='black')
    ax.add_patch(box)
    if text:
        ax.text(x, y, text, ha='center', va='center', fontsize=8)

# Function to draw arrows
def draw_arrow(ax, start, end, color='black'):
    arrow = FancyArrowPatch(start, end, arrowstyle='->', mutation_scale=15, 
                            color=color)
    ax.add_patch(arrow)

# Create figure
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

# Draw Custom GNN container box
draw_box(ax, 0.5, 0.5, 0.5, 0.8, 'Custom GNN', color='none', fill=False)

# GNN subcomponents
sub_labels = ['Encoder', 'GCN Layers', 'Attention', 'Expression Projector', 
              'Differentiable GPR', 'FBA Layer']
sub_ys = [0.85, 0.75, 0.65, 0.55, 0.45, 0.35]
for label, y in zip(sub_labels, sub_ys):
    draw_box(ax, 0.5, y, 0.2, 0.05, label)

# Arrows between GNN subcomponents
for i in range(len(sub_labels)-1):
    start_y = sub_ys[i] - 0.025
    end_y = sub_ys[i+1] + 0.025
    draw_arrow(ax, (0.5, start_y), (0.5, end_y))

# Input boxes
draw_box(ax, 0.1, 0.6, 0.15, 0.05, 'node_features')
draw_box(ax, 0.1, 0.5, 0.15, 0.05, 'edge_index')
draw_box(ax, 0.1, 0.4, 0.15, 0.05, 'distances')
draw_box(ax, 0.1, 0.2, 0.15, 0.05, 'Metabolic Model')

# Arrows from inputs to GNN components
draw_arrow(ax, (0.175, 0.6), (0.25, 0.85))  # node_features to Encoder
draw_arrow(ax, (0.175, 0.5), (0.25, 0.75))  # edge_index to GCN Layers
draw_arrow(ax, (0.175, 0.4), (0.25, 0.65))  # distances to Attention
draw_arrow(ax, (0.175, 0.2), (0.25, 0.45))  # Metabolic Model to Differentiable GPR

# Output and loss
draw_box(ax, 0.8, 0.5, 0.15, 0.05, 'Predicted Fluxes')
draw_box(ax, 0.8, 0.6, 0.15, 0.05, 'Experimental Fluxes')
draw_box(ax, 0.95, 0.5, 0.15, 0.05, 'Supervised Loss')

# Arrows from FBA Layer to Predicted Fluxes
draw_arrow(ax, (0.5, 0.35 - 0.025), (0.8 - 0.075, 0.5))

# Arrows to Supervised Loss
draw_arrow(ax, (0.8 + 0.075, 0.5), (0.95 - 0.075, 0.5))  # Predicted Fluxes
draw_arrow(ax, (0.8 + 0.075, 0.6), (0.95 - 0.075, 0.5))  # Experimental Fluxes

# Training note
ax.text(0.5, 0.05, 'Trained using supervised loss', ha='center', fontsize=10)

# Save the diagram
plt.savefig('model_architecture.png')
plt.close()

print("Schematic diagram saved as 'model_architecture.png'")