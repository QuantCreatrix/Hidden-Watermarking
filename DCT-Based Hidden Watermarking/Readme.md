# DCT-Based Image Watermarking (MATLAB)

## Overview

This project implements a **basic DCT-based image watermarking algorithm** using:

* Block-wise 8×8 DCT
* Mid-frequency coefficient embedding
* Pseudorandom permutation for security

The watermark is embedded into the host image and later extracted to evaluate similarity using **Normalized Correlation (NC)**.

---

## Features

* Grayscale image watermarking
* Mid-frequency DCT embedding
* Key-based watermark scrambling
* Watermark extraction and similarity measurement

---

## Requirements

* MATLAB (R2018 or later recommended)
* Image Processing Toolbox

---

## How to Run

1. Place the following files in the same folder:

   * Main MATLAB script
   * `kohli.jpeg` (input image)
   * `Watermark.png` (binary watermark)

2. Open MATLAB and navigate to the folder

3. Run the script:

```matlab id="run1"
main_script_name
```

---

## Output

* Original image
* Watermarked image
* Extracted watermark
* Normalized Correlation (NC) value

---

## Notes

* Watermark size is automatically resized to `(M/8 × N/8)`
* Embedding is done in the **(4,4) DCT coefficient** of each block
* Embedding strength is fixed (`±0.05`)
* Random permutation key ensures basic security

---
