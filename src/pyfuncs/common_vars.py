

BASE_DIR = "/mnt/data/Proyectos/kranocito"
SEED = 0

CELLBENDER_FIXED_ARGS = [
    "--epochs", "150",
    "--fpr", "0.05",     
    "--learning-rate", "1e-4",
    "--low-count-threshold", "5",
    "--z-dim", "64",
    "--z-layers", "512",
    "--empty-drop-training-fraction", "0.2",
    "--projected-ambient-count-threshold", "0.1",
    "--checkpoint-mins", "15",
    "--cuda",
]