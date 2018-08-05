clear; clc; close all;
%%
bs=4;
im = input('Enter path to image: ');
out = [im(1:end-4) '_compressed.' im(end-2:end)];
btcode2(im, bs, bs, out);
