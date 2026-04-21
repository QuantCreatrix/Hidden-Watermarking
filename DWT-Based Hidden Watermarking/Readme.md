\# 🧠 Blind Medical Image Watermarking (LBP-DWT)



\## 📌 Overview



This project implements a \*\*blind watermarking technique\*\* for medical images based on:



\* \*\*DWT (Discrete Wavelet Transform)\*\*

\* \*\*LBP (Local Binary Pattern)\*\*

\* \*\*Arnold Transform (for scrambling)\*\*

\* \*\*Redundant embedding with majority voting\*\*



The method ensures:



\* High \*\*imperceptibility\*\* (measured using PSNR, SSIM)

\* Strong \*\*robustness\*\* against common attacks (noise, compression, filtering, cropping)



\---



\## 📂 File



\* Main MATLAB script:





\---



\## ⚙️ Requirements



\* MATLAB (R2020 or later recommended)

\* Toolboxes:



&#x20; \* Image Processing Toolbox

&#x20; \* (Optional) Wavelet Toolbox



\---



\## ▶️ How to Run



\### Step 1: Prepare Image



\* Place your medical image in the same folder

\* Rename it as:



```

medical1.png

```



(Or modify `img\_path` in the code)



\---



\### Step 2: Run the Script



In MATLAB:



```matlab

run('your\_script\_name.m')

```



\---



\### Step 3: What Happens



The script will automatically:



1\. Load the image (or generate a phantom if not found)

2\. Create a binary watermark

3\. Embed watermark using LBP-DWT method

4\. Apply multiple attacks

5\. Extract watermark from attacked images

6\. Display:



&#x20;  \* Original vs Watermarked image

&#x20;  \* Extracted watermarks

&#x20;  \* NC \& BCR graphs



\---



\## 📊 Output Metrics



\* \*\*PSNR\*\* → Image quality after embedding

\* \*\*SSIM\*\* → Structural similarity

\* \*\*NC (Normalized Correlation)\*\* → Watermark similarity

\* \*\*BCR (Bit Correct Rate)\*\* → Bit accuracy



\---



\## 🔧 Key Parameters



You can tune these in the script:



```matlab

ARNOLD\_KEY = 5;   % Security (scrambling)

REDUNDANCY = 3;   % Robustness (higher = more reliable)

DWT\_LEVEL  = 2;   % Wavelet depth

```



\---



\## ⚠️ Notes



\* Watermark is generated inside the script (you can replace it)

\* Image size is automatically resized to \*\*512×512\*\*

\* Extraction is \*\*blind\*\* (original image not required)



\---



\## 🚀 Possible Improvements



\* Add geometric attack handling (rotation/scale recovery)

\* Use adaptive embedding strength

\* Replace synthetic watermark with real logo/text



\---





