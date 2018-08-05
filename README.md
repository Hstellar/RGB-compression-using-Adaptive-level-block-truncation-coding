#RGB image compression using Adaptive level Block truncation coding
Multi-level: 2 or 4 level Block Truncation Coding: A lossy compression technique for image compression Adaptive compression: each block is compressed according to the level selected

##Algorithm for level selection
1.Read image 2.Extract RGB values of each pixel 3.Find MAX and MIN values of R plane 4.Compute threshold as average of max and min values 5.Divide R matrix into blocks of 4x4 6.Find MAX and MIN values of each blocks 7.Select 4 level BTC if MAX-MIN >= threshold else select 2 level BTC 8.Similarly perform steps 3-7 for each plane

##Algorithm (2 Level Block Truncation Coding)
1.Calculate mean of each 4x4 block 2.Generate bitmap based on the means 3.Calculate RV0, RV1 (RV: representative value)(formula given in code) 4.Construct new (compressed) matrix based on RV0, RV1 values

##Algorithm (4 Level Block Truncation Coding)
1.Calculate M, R (dynamic range values for threshold calculation) for each 4x4 block(formula given in code) 2.Calculate 3 threshold values, T1, T2, T3, for each RGB planes(formula given in code) 3.Generate bitmap based on T1, T2, T3 values i.e if pixel value > T3 then map it to 3, if T2<= pixel value >= T3 then map it to 2 , if T1<= pixel value >= T2 then map it to 1, else 0.
