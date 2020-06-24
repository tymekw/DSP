clear;
close all;
img = double(imread('lena512.png'));
q=80;
bits = jpegCode(img,q);