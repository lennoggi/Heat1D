# Plotting utilities
1. Generate evolution frames:
   ```
   python3 Plot.py3
   ```
2. Edit the `paths` list in `ffmpeg_list_generator.sh`
3. Generate the list of frames with their time duration for `ffmpeg` (see step 4.):
   ```
   ./ffmpeg_list_generator.sh
   ```
4. Generate the movie:
   ```
   ffmpeg -f concat -safe 0 -i ffmpeg_list.txt -c:v libx264 <outname>.mp4
   ```
