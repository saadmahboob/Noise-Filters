## Noise Reduction Filters

This floder contains noise reduction algorithim. The recorded noisy data for real situations is also avaiable.  

One Matlab function `mimranh.m` is noise reduction alg based on Spectral Subtraction based on Boll-79  
It includes  
1. Amplitude Spectral Subtraction  
2. Magnitude Averaging  
3. Residual Noise Reduction  

Second Matlab function `Wiener_KIST.m` is based on Wiener Noise Suppressor with TSNR & HRNR algorithms  
This Wiener filter is based on tracking a priori SNR using Decision-Directed method, proposed by Plapous et. al. 2006. The two-step noise reduction (TSNR) technique removes the annoying reverberation effect while maintaining the benefits of the decision-directed approach. However, classic short-time noise reduction techniques, including TSNR, introduce harmonic distortion in the enhanced speech. To overcome this problem, a method called harmonic regeneration noise reduction (HRNR)is implemented in order to refine the a priori SNR used to compute a spectral gain able to preserve the speech harmonics.  

Reference: Plapous, C.; Marro, C.; Scalart, P., "Improved Signal-to-Noise Ratio Estimation for Speech Enhancement", IEEE Transactions on Audio, Speech, and Language Processing, Vol. 14, Issue 6, pp. 2098-2108, Nov.2006  

For more detail please visit documentation page by ![CLICK HERE](https://github.com/mimranh/Noise-Filters/wiki)
