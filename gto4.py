import moviepy.editor as mp

clip = mp.VideoFileClip("./mp4/Sq_stable_n_1_alpha_0.001.gif")
clip.write_videofile("./mp4/sq_1_0.001.mp4")