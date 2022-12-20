import moviepy.editor as mp

clip = mp.VideoFileClip("./gto4/sq_1_001_long_5000.gif")
clip.write_videofile("./gto4/sq_1_001_long_5000.mp4")