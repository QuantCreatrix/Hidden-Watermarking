# LGAD-W: LBP-Guided Adaptive DCT Watermarking

## Overview

This project implements a **robust image watermarking technique** using:

* DCT (Discrete Cosine Transform)
* LBP-based texture analysis
* Adaptive block selection (variance + texture)
* Differential embedding strategy

The method ensures **high PSNR, good SSIM, and strong robustness** against attacks like JPEG compression, noise, filtering, and cropping.

---

## Features

* Adaptive watermark embedding using two-factor block selection
* Arnold transform for watermark scrambling
* Robustness testing against multiple attacks
* Performance metrics: **PSNR, SSIM, NC, BCR**
* Automatic result visualization and saving

---

## Requirements

* MATLAB (R2018 or later recommended)
* Image Processing Toolbox

---

## How to Run

1. Place the main script and all functions in the same folder
2. Add an input image (default: `cameraman.tif`)
3. Open MATLAB and navigate to the project folder
4. Run the script:

```matlab
main_script_name
```

---

## Output

* Watermarked image
* Extracted watermark
* Performance metrics (PSNR, SSIM, NC, BCR)
* Saved figures in `LGADW_Results/`:

  * Overview
  * Block maps
  * Robustness graphs
  * Attack results

---

## Notes

* Watermark size is automatically determined based on image capacity
* Parameters like `ALPHA`, `VARIANCE_PERCENTILE`, and `LBP_PERCENTILE` can be tuned for performance
* Designed for grayscale images (RGB converted automatically)

---

