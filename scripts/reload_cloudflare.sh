conda activate single-cell

tmux new-session

# En la nueva sesion
cellxgene launch /mnt/data/Proyectos/kranocito/data/ARAUZO_03/processed_adatas/combined_matrix_cellxgene_v2.h5ad --host 127.0.0.1 --port 5005 --disable-annotations --disable-diffexp

# Ctrl + b + c
cloudflared tunnel --protocol http2 --url http://127.0.0.1:5005 --loglevel info 



# repeat for oprescu
cellxgene launch /mnt/data/Proyectos/kranocito/data/processed/oprescu_processed_clean_FAP_with_phate.h5ad --host 127.0.0.1 --port 5006 --disable-annotations --disable-diffexp
cloudflared tunnel --protocol http2 --url http://127.0.0.1:5006 --loglevel info 

 
