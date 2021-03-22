close all;
clear;

[sIn, Fs] = audioread("gas_station.wav");

# 1 s -> 16000
# 0.010 s -> 160
# 0.025 s -> 400

FRAME_SIZE = 400;
SHIFT_SIZE = 160;

frameStart = [];

for i = 1:SHIFT_SIZE:length(sIn)
  frameStart(length(frameStart)+1) = i;
end

# Pad
for i = length(sIn)+1:frameStart(length(frameStart))+FRAME_SIZE
  sIn(i) = 0;
end

# Calculate sOf
sOf = zeros(length(sIn), 1);
sOf(1) = sIn(1);
for i = 2:length(sIn)
  sOf(i) = sIn(i) - sIn(i-1) + 0.999 * sOf(i-1);
end

# Calculate log frame energy
logE = [];
for fi = 1:length(frameStart)
  currsum = 0;
  for i = frameStart(fi):frameStart(fi)+FRAME_SIZE-1
    currsum += sOf(i) ** 2;
  endfor
  logE(length(logE)+1) = max(-50, log(currsum));
  
  if fi == 50
    logEframe50 = max(-50, log(currsum));
  endif
end

# Calculate sPe
sPe = zeros(length(sIn), 1);
sPe(1) = sOf(1);
for i = 2:length(sIn)
  sPe(i) = sOf(i) - 0.97 * sOf(i-1);
end

# Plot sPe (magnitude)
# pre-emphasis is high pass filter because it boost the magnitude of high frequency and make low frequency magnitude much lower
plot(abs(fft(sOf))); # Before
plot(abs(fft(sPe))); # After

# Extract frames
NDFT = 512;
framesOf = zeros(length(frameStart), FRAME_SIZE);
framesPe = zeros(length(frameStart), FRAME_SIZE);
framesHam = zeros(length(frameStart), FRAME_SIZE);
framesFFT = zeros(length(frameStart), NDFT);
for fi = 1:length(frameStart)
  for i = frameStart(fi):frameStart(fi)+FRAME_SIZE-1
    framesOf(fi, i-frameStart(fi)+1) = sOf(i);
    framesPe(fi, i-frameStart(fi)+1) = sPe(i);
    aa = framesPe(fi, :);
    aa_w = framesHam(fi, :) = hamming(length(aa))'.*aa;
    framesFFT(fi, :) = fft(aa_w,NDFT);
  endfor
endfor

# Verify framing correct
currsum = 0
for i = 1:FRAME_SIZE
  currsum += framesOf(50, i) ** 2;
endfor
logEframe50New = max(-50, log(currsum));

# Plot frame 50 fft
f = 0:Fs/NDFT:Fs-Fs/NDFT;
figure;
hold on;
AA = framesFFT(50, :);
plot(f,abs(AA)/max(abs(AA)));
axis([0,8000,0,1]);
grid on;
xlabel("f (Hz.)");

# Load mel
load("mel_filters.mat");
plot((0:256)'/256*8000, mel_filters);
title("Mel filter banks");
xlabel("frequency (Hz)");
ylabel("weight");
