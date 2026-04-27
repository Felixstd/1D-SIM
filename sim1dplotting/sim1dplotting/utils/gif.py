import glob
from PIL import Image
import numpy as np
import cv2
import re

"""
Code used to make gifs
"""

def make_gif(frame_folder):
    frames = [Image.open(image, 'r') for image in sorted(glob.glob(frame_folder+"nonviscid_wave_*"))]#, #key = lambda x: int(x.replace('age_grid_', '')))]
    for image in sorted(glob.glob(frame_folder+"nonviscid_wave_*")) : 
        print(image)
    frame_one = frames[0]
    frame_one.save("inviscid.gif", format="GIF", append_images=frames,
               save_all=True, duration=300, loop=0)
    
make_gif('/aos/home/fstdenis/1D-SIM/nonviscid/')
