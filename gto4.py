import moviepy.editor as mp

clip = mp.VideoFileClip("./gto4/tri_2_N_8000.gif")
clip.write_videofile("./gto4/tri_2_N_8000.mp4")