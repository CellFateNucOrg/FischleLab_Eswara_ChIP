import os
import glob
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from pdf2image import convert_from_path

# Set the root directory
root_dir = '/mnt/external.data/MeisterLab/FischleLab_KarthikEswara/ChIP'

#samples = pd.read_csv(os.path.join(root_dir, 'sampleSheet.csv'))

out_dir="gathered_qc_plots"
os.makedirs(out_dir,exist_ok=True)

########
# deeptools fingerprint plots
sub_dir = 'bwa/merged_library/deepTools/plotFingerprint'
images_info=glob.glob(os.path.join(root_dir, sub_dir, '*plotFingerprint.pdf'))
images_info.sort()  # Ensure consistent order

pngs = []
for pdf in images_info:
    pages = convert_from_path(pdf, dpi=300)
    pngs.append(pages[0])   # assuming 1-page PDFs

# Paginate in chunks of 12
chunks = [pngs[i:i+12] for i in range(0, len(pngs), 12)]


with PdfPages(os.path.join(out_dir, 'plotFingerprint_all.pdf')) as pdf:
    for chunk in chunks:
        fig, axs = plt.subplots(3, 4, figsize=(12, 8))
        axs = axs.T.flatten()
        
        for ax, img in zip(axs, chunk):
            ax.imshow(img)
            ax.axis('off')
        
        for ax in axs[len(chunk):]:
            ax.axis('off')
        
        fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.02, hspace=0.02)
        #plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


##### 
# deeptools heatmaps
sub_dir = 'bwa/merged_library/deepTools/plotProfile'
images_info=glob.glob(os.path.join(root_dir, sub_dir, '*plotHeatmap.pdf'))
images_info.sort()  # Ensure consistent order

pngs = []
for pdf in images_info:
    pages = convert_from_path(pdf, dpi=300)
    pngs.append(pages[0])   # assuming 1-page PDFs

# Paginate in chunks of 18
chunks = [pngs[i:i+12] for i in range(0, len(pngs), 12)]


with PdfPages(os.path.join(out_dir, 'plotHeatmap_all.pdf')) as pdf:
    for chunk in chunks:
        fig, axs = plt.subplots(1, 12, figsize=(12, 8))
        axs = axs.flatten()
        
        for ax, img in zip(axs, chunk):
            ax.imshow(img)
            ax.axis('off')
        
        for ax in axs[len(chunk):]:
            ax.axis('off')
        
        fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.02, hspace=0.02)
        #plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

######
# phantompeakqualtools
sub_dir = 'bwa/merged_library/phantompeakqualtools'
images_info=glob.glob(os.path.join(root_dir, sub_dir, '*.spp.pdf'))
images_info.sort()  # Ensure consistent order

pngs = []
for pdf in images_info:
    pages = convert_from_path(pdf, dpi=300)
    pngs.append(pages[0])   # assuming 1-page PDFs

# Paginate in chunks of 18
chunks = [pngs[i:i+12] for i in range(0, len(pngs), 12)]


with PdfPages(os.path.join(out_dir, 'phantompeakqualtools_all.pdf')) as pdf:
    for chunk in chunks:
        fig, axs = plt.subplots(3, 4, figsize=(12, 8))
        axs = axs.T.flatten()
        
        for ax, img in zip(axs, chunk):
            ax.imshow(img)
            ax.axis('off')
        
        for ax in axs[len(chunk):]:
            ax.axis('off')
        
        fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.02, hspace=0.02)
        #plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)