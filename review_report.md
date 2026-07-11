# Báo cáo Review lần 2: Bài tập lớn Phương pháp số cho học máy

> **Đề tài**: Các phương pháp tối ưu chuyên biệt trong học máy (Chapter 12.1, 12.2, 12.3)
> **Nhóm**: 5 — Đại học Phenikaa
> **Tài liệu tham chiếu**: *Numerical Algorithms* — Justin Solomon (2015), Chapter 12

---

## 1. Các cải tiến đã thực hiện so với lần review 1

| Vấn đề từ lần 1 | Trạng thái | Ghi chú |
|---|---|---|
| Lỗi `\documentclass` khai báo 2 lần | ✅ Đã sửa | Dòng 48 main.tex đã gọn |
| 2 trang bìa trùng lặp | ✅ Đã sửa | Gộp thành 1 trang bìa duy nhất |
| Tên sách tham chiếu không rõ | ✅ Đã sửa | Đổi sang cite Solomon2015 với tên sách đầy đủ |
| Thiếu chi tiết ADMM convergence | ✅ Đã bổ sung | Thêm subsubsection "Nhận xét về hội tụ và tham số ρ" (sec4, dòng 233–240) |
| Convention dấu Lagrangian khác tài liệu | ✅ Đã bổ sung | Thêm "Chú ý" giải thích convention (sec4, dòng 176–178) |
| Phần thực nghiệm rỗng | ⚠️ Cải thiện một phần | Có đề cập notebook, nhưng **chưa có hình/bảng chèn vào báo cáo** |
| Không có hình minh họa nào | 🔴 Chưa sửa | Vẫn không có figure nào trong báo cáo |
| Section 5 trùng lặp nội dung | 🟡 Chưa sửa | Vẫn lặp lại công thức từ section 2–4 |
| Thiếu MSSV 2 thành viên | 🔴 Chưa sửa | Dòng 131–132 vẫn để trống |

---

## 2. Đánh giá Coverage (mức độ bao phủ kiến thức)

### 2.1. Đối chiếu từng đầu mục với tài liệu gốc

| Nội dung trong Chapter 12 (Solomon 2015) | Báo cáo | Chi tiết |
|---|:---:|---|
| **12.1 Nonlinear Least-Squares** | ✅ | Section 2, đầy đủ |
| — Ví dụ exponential regression y = ce^ax | ✅ | Sec 2.1 |
| — Công thức NLS tổng quát, Jacobian, Hessian | ✅ | Sec 2.2 |
| — 12.1.1 Gauss-Newton: linearization, normal eq, update | ✅ | Sec 2.3 |
| — GN convergence: quadratic near optimum, fails when E* large | ✅ | Sec 2.3 |
| — 12.1.2 Levenberg-Marquardt: trust region, KKT, damping | ✅ | Sec 2.4 |
| — LM ↔ GN khi λ nhỏ, LM ↔ GD khi λ lớn | ✅ | Sec 2.4 |
| — Adaptive λ adjustment | ✅ | Sec 2.4 (bổ sung tỷ số ρ_k) |
| **12.2 IRLS** | ✅ | Section 3, đầy đủ |
| — Công thức tổng quát E = Σ f_i [g_i]² | ✅ | Sec 3.1 |
| — Ex 12.1: Lp optimization | ✅ | Sec 3.2 |
| — Ex 12.2: L1 optimization + δ smoothing | ✅ | Sec 3.3 |
| — Ex 12.3: Weiszfeld (geometric median) | ✅ | Sec 3.4 |
| — Nhận xét hội tụ IRLS, L1 convergence | ✅ | Sec 3.5 |
| **12.3 Coordinate Descent and Alternation** | ✅ | Section 4, đầy đủ |
| — Monotone decrease property | ✅ | Sec 4.1 |
| — 12.3.1 Identifying candidates for alternation | ✅ | Sec 4.2 |
| — Ex 12.4: Generalized PCA / ALS | ✅ | Sec 4.3 |
| — Ex 12.5: ARAP (Procrustes + sparse linear) | ✅ | Sec 4.3 |
| — Ex 12.6: Coordinate descent cho LS | ✅ | Sec 4.3 |
| — Ex 12.7: k-means clustering | ✅ | Sec 4.3 |
| — 12.3.2 Augmented Lagrangians and ADMM | ✅ | Sec 4.4 |
| — Quadratic penalty → Lagrangian → Aug. Lagrangian | ✅ | Sec 4.4 |
| — ADMM general form: 3-step update | ✅ | Sec 4.4 |
| — ADMM convergence guarantees (convex case) | ✅ | **Mới bổ sung** |
| — ADMM: fast approximate, slow to high accuracy | ✅ | **Mới bổ sung** |
| — ρ selection strategy | ✅ | **Mới bổ sung** |
| — Ex 12.8: NNLS via ADMM | ✅ | Sec 4.4 |
| — Ex 12.9: Geometric median via ADMM | ✅ | Sec 4.4 |

> [!TIP]
> **Coverage: 10/10** — Tất cả nội dung trong 12.1, 12.2, 12.3 đều được trình bày. Phần ADMM convergence đã được bổ sung đầy đủ.

---

## 3. Đánh giá Correctness (tính chính xác)

### 3.1. Kiểm tra công thức — kết quả

| Công thức | Vị trí | Kết quả |
|---|---|:---:|
| E_NLS = 1/2 Σ [f_i]² | sec2 L32-33 | ✅ Khớp |
| ∇E = J^T F | sec2 L58 | ✅ Khớp |
| Hessian = J^TJ + Σ f_i ∇²f_i | sec2 L63-64 | ✅ Khớp |
| GN linearization F(x+δ) ≈ F + Jδ | sec2 L81 | ✅ Khớp |
| GN normal eq: J^TJδ = -J^TF | sec2 L103 | ✅ Khớp |
| GN update: x_{k+1} = x_k - (J^TJ)^{-1}J^TF | sec2 L113-114 | ✅ Khớp |
| LM trust region formulation | sec2 L141-151 | ✅ Khớp |
| KKT → (H + λI)δ = -g | sec2 L176-181 | ✅ Khớp |
| LM: (J^TJ + λI)δ = -J^TF | sec2 L186-187 | ✅ Khớp |
| LM ↔ GD: δ ≈ -(1/λ)∇E | sec2 L208-209 | ✅ Khớp |
| IRLS: E = Σ f_i [g_i]² | sec3 L11-13 | ✅ Khớp |
| Lp: |r|^p = |r|^{p-2} r² | sec3 L65-67 | ✅ Khớp |
| L1 weight: 1/max(|r|, δ) | sec3 L106-108 | ✅ Khớp |
| Weiszfeld: x = Σw_i x_i / Σw_i | sec3 L174-177 | ✅ Khớp |
| Alternation monotone decrease | sec4 L18-27 | ✅ Khớp |
| PCA as ALS: min ‖X-CY‖² | sec4 L48 | ✅ Khớp |
| ARAP energy + rotation constraint | sec4 L69-77 | ✅ Khớp |
| CD update formula | sec4 L110-113 | ✅ Khớp |
| k-means objective + 2-step | sec4 L120-147 | ✅ Khớp |
| Augmented Lagrangian L_ρ | sec4 L183-185 | ✅ Convention +λ (có ghi chú) |
| ADMM 3-step | sec4 L217-228 | ✅ Khớp |
| NNLS ADMM: x, z, λ updates | sec4 L278-303 | ✅ Khớp |
| Geo-median ADMM: x, z, λ updates | sec4 L322-348 | ✅ Khớp |
| Lasso ADMM: soft-thresholding | sec5 L114-141 | ✅ Khớp |

### 3.2. Hai điểm khác biệt (không phải lỗi)

> [!NOTE]
> **Điểm 1: Lp weight smoothing**
> Sách gốc viết trực tiếp f_i = (a_i·x - b_i)^{p-2}. Báo cáo dùng phiên bản stabilized w_i = ((r_i)² + ε²)^{p/2-1}. Đây là cải tiến đúng cho triển khai số. Nên thêm 1 câu ghi chú so sánh.

> [!NOTE]
> **Điểm 2: Tỷ số ρ_k trong LM**
> Không có trong sách Solomon — lấy từ literature trust-region. Nội dung đúng. Nên thêm \cite{NocedalWright2006} ở chỗ này.

### 3.3. Kết luận

> [!TIP]
> **Không phát hiện lỗi sai công thức nào.** Convention Lagrangian đã được giải thích rõ ràng. Tất cả công thức khớp với tài liệu gốc.

---

## 4. Đánh giá diễn đạt và trình bày

### 4.1. Điểm mạnh

- Cấu trúc logic rõ ràng, flow tốt
- Giải thích trực giác vượt mức liệt kê công thức
- Convention sign Lagrangian được ghi chú chuyên nghiệp
- ADMM convergence + ρ strategy đầy đủ
- BibTeX sạch, tài liệu tham khảo chính xác

### 4.2. Vấn đề còn tồn đọng

> [!CAUTION]
> **Vấn đề 1 (nghiêm trọng nhất): Thực nghiệm không có kết quả trực quan**
>
> Section 6 đề cập notebook nhưng **không chèn hình, bảng hay số liệu nào**. Nếu giáo viên mở PDF và không thấy kết quả thực nghiệm, phần này gần như không có giá trị điểm.

> [!IMPORTANT]
> **Vấn đề 2: Không có figure nào (ngoài logo trường)**
>
> ~30 trang báo cáo toán học mà không hình minh họa. Tài liệu gốc có Figure 12.1 (CD zigzag), 12.2 (k-means), 12.3 (Aug. Lagrangian).

> [!NOTE]
> **Vấn đề 3: Section 5 lặp công thức** — nên dùng \eqref trở lại section gốc.

> [!NOTE]
> **Vấn đề 4: MSSV thiếu 2 người** — main.tex dòng 131-132.

> [!NOTE]
> **Vấn đề 5: `\usepackage{amsmath}` load 2 lần** — dòng 8 và dòng 93.

---

## 5. Đánh giá theo tiêu chí (so sánh lần 1 vs lần 2)

| Tiêu chí | Lần 1 | Lần 2 | Nhận xét |
|-----------|:---:|:---:|---------|
| **Bao phủ đầu mục** | 9.5 | **10** | Đã bổ sung ADMM convergence, ρ strategy |
| **Chính xác công thức** | 9.5 | **9.5** | Vẫn chính xác, không lỗi mới |
| **Diễn đạt / giải thích** | 8.5 | **9** | Convention note, ADMM nhận xét đầy đủ |
| **Thực nghiệm** | 5 | **5.5** | Có đề cập notebook nhưng không chèn hình/bảng |
| **Hình thức trình bày** | 8 | **8.5** | Sửa documentclass, gộp bìa, bib chuẩn |

---

## 6. Đánh giá tổng thể

### Ước tính điểm hiện tại: 8 – 8.5/10

Cải thiện rõ rệt so với lần 1 (từ ~7.5-8 lên ~8-8.5). Phần lý thuyết đã rất tốt (9.5-10/10). Điểm trừ chính vẫn là phần thực nghiệm thiếu kết quả và không có hình minh họa.

### Để đạt 9/10, cần làm thêm (ưu tiên giảm dần):

**1. Chèn hình thực nghiệm vào báo cáo** (ưu tiên cao nhất)

Chạy notebook, export hình, thêm vào sec6 với `\includegraphics`. Cần ít nhất 3-4 hình:
- GN/LM fit + convergence plot
- IRLS L1 vs L2 khi có ngoại lai
- k-means clustering result
- ADMM primal residual convergence

**2. Thêm hình minh họa lý thuyết**

Tối thiểu 2 hình:
- Coordinate descent zigzag (tương tự Figure 12.1 sách)
- Augmented Lagrangian penalty levels (tương tự Figure 12.3)

**3. Rút gọn Section 5** — dùng \eqref thay vì lặp công thức

**4. Điền MSSV** cho 2 thành viên còn thiếu

**5. Sửa lỗi nhỏ** — xóa `\usepackage{amsmath}` dòng 93, thêm cite cho ρ_k

---

## 7. Bảng tổng hợp vấn đề còn lại

| # | Mức độ | Vấn đề | Hành động |
|---|--------|--------|-----------|
| 1 | 🔴 Nghiêm trọng | Thực nghiệm không có hình/bảng trong PDF | Chạy notebook → chèn hình |
| 2 | 🟠 Quan trọng | Không có figure nào (ngoài logo) | Thêm ≥ 5 hình |
| 3 | 🟡 Trung bình | Section 5 lặp công thức | Dùng \eqref |
| 4 | 🔵 Nhẹ | MSSV thiếu 2 người | Điền MSSV |
| 5 | 🔵 Nhẹ | amsmath load 2 lần | Xóa dòng 93 |
| 6 | 🔵 Nhẹ | Thiếu cite cho ρ_k ratio | Thêm cite Nocedal |

---

*Review lần 2 — đối chiếu với Chapter 12 (12.1–12.3), Numerical Algorithms, Justin Solomon (2015)*
