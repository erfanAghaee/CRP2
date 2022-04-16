ffmpeg -framerate 8 -pattern_type glob -i '*.png' \
  -c:v libx264 -pix_fmt yuv420p out.mp4

# ffmpeg -framerate 5 -i nets_%d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4