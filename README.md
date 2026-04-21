# Hidden Watermarking: Research, Implementation & Proposed Method

## Overview

This repository presents a **structured study and implementation of digital image watermarking techniques**, culminating in a **novel proposed method (LGAD-W)**.

The work progresses through:

1. Classical DCT-based watermarking
2. Modern hybrid LBP–DWT watermarking
3. Proposed **LGAD-W (LBP-Guided Adaptive DCT)** method

The objective is to analyze, compare, and improve the **trade-off between imperceptibility and robustness**.

---

## 1. Classical Approach – DCT Watermarking

Based on early work in digital watermarking 

### Key Ideas:

* Block-based **8×8 DCT transformation**
* Watermark embedded in **mid-frequency coefficients**
* Use of **pseudorandom permutation** for security

### Strengths:

* Good compatibility with JPEG compression
* Simple and computationally efficient

### Limitations:

* Limited robustness to strong attacks
* Fixed embedding strategy

---

## 2. Modern Approach – LBP–DWT Watermarking

Based on recent research for secure image transmission 

### Key Ideas:

* **DWT decomposition (LL, LH, HL, HH sub-bands)**
* **LBP (Local Binary Pattern)** for texture analysis
* Embedding using **LSB in LL sub-band**
* **Arnold Transform** for watermark scrambling

### Strengths:

* Better robustness against noise, compression, filtering
* Maintains high image quality (high PSNR & SSIM)

### Limitations:

* Higher computational complexity
* LSB-based embedding may still be vulnerable

---

## 3. Proposed Method – LGAD-W

Developed in this project 

### Core Contributions:

* **Two-factor block selection**

  * Block variance (energy)
  * LBP texture complexity
* **Differential DCT embedding**

  * Encodes bits using coefficient differences
  * Reduces distortion and improves robustness
* **Arnold scrambling**

  * Enhances security and error distribution

---

### Why LGAD-W is Better

| Aspect           | Classical DCT | LBP–DWT  | LGAD-W (Proposed)      |
| ---------------- | ------------- | -------- | ---------------------- |
| Robustness       | Medium        | High     | **Very High**          |
| Imperceptibility | Medium        | High     | **High (PSNR > 40dB)** |
| Complexity       | Low           | High     | Moderate               |
| Security         | Basic         | Moderate | **Enhanced**           |

---

## Features Across Implementations

* Watermark embedding and extraction
* Performance metrics:

  * PSNR (Peak Signal-to-Noise Ratio)
  * SSIM (Structural Similarity Index)
  * NC (Normalized Correlation)
  * BCR (Bit Correct Rate)
* Robustness testing:

  * JPEG compression
  * Gaussian noise
  * Salt & Pepper noise
  * Filtering
  * Cropping

---

## How to Run

### Step 1: Open MATLAB

Navigate to any method folder:

```
Paper1_DCT / Paper2_LBP_DWT / Proposed_LGADW
```

### Step 2: Run the main script

```matlab
main_script_name
```

### Step 3: Outputs

* Watermarked image
* Extracted watermark
* Performance metrics
* (Proposed method) robustness analysis plots

---

## Objective

This repository aims to:

* Study the evolution of watermarking techniques
* Compare classical and modern approaches
* Develop a **more robust and adaptive watermarking method**

---

## Key Takeaway

Watermarking is fundamentally a **trade-off problem**:

* Too strong → visible distortion
* Too weak → poor robustness

LGAD-W improves this balance by combining:

* **Frequency domain embedding (DCT)**
* **Texture awareness (LBP)**
* **Adaptive block selection**

---

## Contributors

* Dhananjay Agrahari
* Anubhav Rathore
* Sandesh Jat
* Aman Kanaujiya

---

## Notes

* Each folder is independent and runnable
* Input images may need to be added manually
* Parameters can be tuned for experimentation

---
