# Newton Forward & Backward Interpolation (Qt GUI)

This is a C++ application with Qt GUI and [muParser](https://beltoforion.de/en/muparser/) library that demonstrates **Newton's Forward and Backward Interpolation** methods with:

- Interactive table for data input
- Polynomial equation generation
- Evaluation of the polynomial at a given `x`
- Newton-Raphson method to solve for `x` given `y`
- Visual switching between Forward and Backward methods

---

## 🖼️ Screenshot

![App Screenshot](screenshot.png)  
*Add your screenshot in the repo with name `screenshot.png` or edit this line accordingly.*

---

## 🚀 Features

- Create a data table with `X` and `f(X)`
- Compute divided differences and generate interpolation polynomial
- Calculate `f(x)` at any input `x`
- Calculate `x` from `f(x)` using Newton-Raphson
- Choose between **Forward** and **Backward** interpolation

---

## 🧰 Requirements

- Qt 5 or 6 (tested on Qt 5.15+)
- muParser library

---

## 🛠️ How to Build and Run

1. **Clone the repo:**

```bash
git clone https://github.com/M-Maaly/newton-forward-backward.git
cd newton-forward-backward
