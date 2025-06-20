# DSP 專案集

本倉庫具開發的三個 MATLAB 數位訊號處理專案，涵蓋 FIR 濾波器頻率響應、陷波濾波器去除雜訊，以及 IIR 濾波器組音符偵測。

---

## 目錄

1. [前置需求](#前置需求)  
2. [檔案結構](#檔案結構)  
3. [專案 1：FIR 濾波器頻率響應與取樣影響](#專案-1fir-濾波器頻率響應與取樣影響)  
4. [專案 2：FIR 陷波濾波器去除雜訊](#專案-2fir-陷波濾波器去除雜訊)  
5. [專案 3：IIR 濾波器組音符偵測](#專案-3iir-濾波器組音符偵測)  
6. [執行方式](#執行方式)  

---

## 前置需求

- MATLAB R2024a 以上  
- Signal Processing Toolbox  

---

## 檔案結構

```

/DSP-Projects
├─ DSP_undergraduate/      # 專案 1：FIR 濾波器頻率響應與取樣影響
│   ├─ DSP_undergraduate.m
│   └─ README.md
├─ DSP_graduate/      # 專案 2：FIR 陷波濾波器去除雜訊
│   ├─ DSP_graduate.m
│   └─ README.md
├─ Final_project/       # 專案 3：IIR 濾波器組音符偵測
│   ├─ toneDetect.m
│   ├─ main.m
│   └─ README.md
└─ README.md                     # 本檔


````

---

## 專案 1：FIR 濾波器頻率響應與取樣影響

### 概述

探索在不同採樣率（如 4000 Hz、8000 Hz 等）下，使用固定「數位頻率」參數設計 FIR 濾波器時，其實際截止頻率隨採樣率變化的現象。

### 實作

1. **濾波器設計**  
   - 使用 `fir1` 或自訂係數設計指定數位截止頻率的 FIR 濾波器。  
2. **頻率響應分析**  
   - 以 `freqz` 繪製頻率響應並擷取歸一化截止頻率。  
3. **訊號測試**  
   - 產生 300 Hz 與 1200 Hz 正弦波，並將其相加後過濾。  
   - 使用 FFT 驗證 1200 Hz 成分的抑制效果。  
4. **採樣率比較**  
   - 在不同採樣率下重複上述流程，並對比實際截止頻率。

| 採樣率 \(f_s\) (Hz) | 歸一化截止頻率 (cycles/sample) | 截止頻率 \(f_c\) (Hz) |
|--------------------:|-------------------------------:|----------------------:|
| 4000               | 0.1348                         | 539.1                 |
| 8000               | 0.1348                         | 1078.1                |

### 結果

- 在 4000 Hz 下能有效抑制 1200 Hz，而保留 300 Hz。  
- 提高至 8000 Hz 時，截止頻率上升至約 1078 Hz，導致對 1200 Hz 的抑制效果下降，凸顯採樣率選擇的重要性。

---

## 專案 2：FIR 陷波濾波器去除雜訊

### 概述

設計三組簡單 3 項式 FIR 陷波濾波器，以去除特定三個雜訊頻率（1575 Hz、3150 Hz、4725 Hz），並保持原始音訊主要內容。

### 實作

1. **頻率識別**  
   - 透過 spectrogram 與 `ginput` 初步選取雜訊帶位置。  
   - 使用 FFT 精確量測各主要雜訊頻率。  
2. **陷波濾波器設計**  
   - 計算數位角頻率 \(\omega = 2\pi f / f_s\)。  
   - 係數 \(A = -2\cos(\omega)\)，使濾波器在該頻率產生零點。  
3. **濾波流程**  
   - 依序對輸入音訊套用三個 3 項 FIR 濾波器。

| 雜訊頻率 (Hz) | \(\omega\) (rad/sample) | 係數 \(A\)   |
|-------------:|-----------------------:|-------------:|
| 1575         | 0.286π                 | −1.246       |
| 3150         | 0.571π                 | −0.442       |
| 4725         | 0.857π                 | +1.802       |

### 結果

- Spectrogram 顯示三條雜訊線成功被抑制。  
- 處理後音訊更為乾淨，且原始內容保留良好。

---

## 專案 3：IIR 濾波器組音符偵測

### 概述

實作一組 13 個 IIR 峰值型帶通濾波器，對應 C4 (261.63 Hz) 至 C5 (523.25 Hz) 每個半音，並根據各濾波器輸出能量及 FFT 峰值，自動判別不明音訊的音符。

### 實作

1. **濾波器設計**  
   - 使用 `iirpeak` 為每個音符頻率設計帶通濾波器。  
2. **偵測流程**  
   1. 讀入 WAV 檔至 `toneDetect(xx, fs)`。  
   2. 將訊號同時輸入所有 13 支濾波器。  
   3. 計算各支濾波器輸出能量並正規化。  
   4. 執行 FFT 找出顯著頻率峰值。  
   5. 根據能量分布與熵值，分類為「精準匹配」、「失諧」或「超出範圍」。

### 結果範例

| 檔名         | 判別結果      | 備註                         |
|-------------|--------------|-----------------------------|
| Signal01.wav | 超出範圍     | 所有濾波器能量皆偏低         |
| Signal03.wav | 精準匹配 G4  | 第 7 支濾波器能量最高        |
| Signal08.wav | 精準匹配 A4  | 第 9 支濾波器能量最高        |

---

## 執行方式

1. 將本倉庫複製至本機：  
   ```bash
   git clone https://github.com/leeminwei/Digital-Signal-Processing.git
   cd (DSP_graduate/DSP_undergraduate/Final_project)    
    ```

2. 打開 MATLAB，將工作目錄切至欲執行之專案資料夾。
3. 執行對應主程式，如：

   ```matlab
   >> DSP_undergrduate    % 專案 1
   >> DSP_fft             % 專案 2
   >> main          % 專案 3
   ```
4. 依照螢幕提示查看繪圖結果及聆聽處理後音訊。


