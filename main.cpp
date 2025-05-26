/************************************************************************/


#include <QApplication>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QTableWidget>
#include <QHeaderView>
#include <QInputDialog>
#include <QMessageBox>
#include <QPainter>
#include <QLabel>
#include <QComboBox>
#include <muParser.h>
#include <cmath>
#include <sstream>
#include <iomanip>

class DividedDifferencesApp : public QWidget {
    Q_OBJECT

public:
    DividedDifferencesApp(QWidget *parent = nullptr) : QWidget(parent) {
        setupUI();
    }

protected:
    void paintEvent(QPaintEvent *event) override {
        QPainter painter(this);
        QLinearGradient gradient(0, 0, width(), height());
        gradient.setColorAt(0, QColor(255, 255, 255));
        gradient.setColorAt(1, QColor(200, 230, 255));
        painter.fillRect(rect(), gradient);
        QWidget::paintEvent(event);
    }

private:
    QTableWidget *tableWidget;
    QPushButton *createButton, *computeButton, *calcButton, *calcXButton, *clearButton;
    QComboBox *methodCombo;
    QLabel *equationLabel;

    void setupUI() {
        QVBoxLayout *layout = new QVBoxLayout(this);

        tableWidget = new QTableWidget(this);
        tableWidget->setEditTriggers(QAbstractItemView::AllEditTriggers);
        tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
        tableWidget->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);
        tableWidget->setStyleSheet("QTableWidget { background-color: #ffffff; }");

        createButton = new QPushButton("Create Table", this);
        computeButton = new QPushButton("Compute Divided Differences", this);
        calcButton = new QPushButton("Calculate Value at X", this);
        calcXButton = new QPushButton("Calculate X from Y", this);
        clearButton = new QPushButton("Clear Table", this);

        QString buttonStyle = "QPushButton { background-color: #6ca0dc; color: white; font-weight: bold; padding: 8px; border-radius: 8px; }"
                              "QPushButton:hover { background-color: #5b8bd0; }";
        createButton->setStyleSheet(buttonStyle);
        computeButton->setStyleSheet(buttonStyle);
        calcButton->setStyleSheet(buttonStyle);
        calcXButton->setStyleSheet(buttonStyle);
        clearButton->setStyleSheet(buttonStyle);

        methodCombo = new QComboBox(this);
        methodCombo->addItem("Forward Interpolation");
        methodCombo->addItem("Backward Interpolation");
        methodCombo->setStyleSheet("QComboBox { padding: 5px; }");

        equationLabel = new QLabel("Polynomial Equation: ", this);
        equationLabel->setStyleSheet("QLabel { font-size: 14px; padding: 10px; }");

        QHBoxLayout *methodLayout = new QHBoxLayout();
        methodLayout->addWidget(new QLabel("Interpolation Method:", this));
        methodLayout->addWidget(methodCombo);

        layout->addWidget(tableWidget);
        layout->addLayout(methodLayout);
        layout->addWidget(equationLabel);
        layout->addWidget(createButton);
        layout->addWidget(computeButton);
        layout->addWidget(calcButton);
        layout->addWidget(calcXButton);
        layout->addWidget(clearButton);

        connect(createButton, &QPushButton::clicked, this, &DividedDifferencesApp::createTable);
        connect(computeButton, &QPushButton::clicked, this, &DividedDifferencesApp::computeDividedDifferences);
        connect(calcButton, &QPushButton::clicked, this, &DividedDifferencesApp::calculateValueAtX);
        connect(calcXButton, &QPushButton::clicked, this, &DividedDifferencesApp::calculateXFromY);
        connect(clearButton, &QPushButton::clicked, this, &DividedDifferencesApp::clearTable);
    }

    double evaluateExpression(const QString &expression) {
        try {
            mu::Parser parser;
            parser.DefineConst(L"pi", M_PI);
            parser.DefineConst(L"e", M_E);
            parser.SetExpr(expression.toStdWString());
            return parser.Eval();
        } catch (mu::Parser::exception_type &e) {
            QMessageBox::warning(this, "Error", QString::fromStdWString(e.GetMsg()));
            return 0.0;
        }
    }

    QString generatePolynomial(const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        int n = x.size();
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);


        ss << "P(x) = ";
        for (int i = 0; i < n; ++i) {
            double coef = f[i][0];
            if (coef == 0) continue;

            if (i > 0) ss << (coef > 0 ? " + " : " - ");
            ss << std::abs(coef);

            for (int j = 0; j < i; ++j) {
                int idx = isForward ? j : (n - 1 - j);
                ss << "(x";
                if (x[idx] >= 0) ss << " - " << x[idx];
                else ss << " + " << -x[idx];
                ss << ")";
            }
        }
        return QString::fromStdString(ss.str());
    }

    double evaluatePolynomial(double XX, const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        int n = x.size();
        double P = 0;
        for (int i = 0; i < n; ++i) {
            double k = 1;
            for (int j = 0; j < i; ++j) {
                int idx = isForward ? j : (n - 1 - j);
                k *= (XX - x[idx]);
            }
            P += k * f[i][0];
        }
        return P;
    }

    double evaluateDerivative(double XX, const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        int n = x.size();
        double deriv = 0;
        for (int i = 1; i < n; ++i) {
            double term = f[i][0];
            for (int j = 0; j < i; ++j) {
                double prod = 1;
                for (int k = 0; k < i; ++k) {
                    if (k == j) continue;
                    int idx = isForward ? k : (n - 1 - k);
                    prod *= (XX - x[idx]);
                }
                deriv += term * prod;
            }
        }
        return deriv;
    }

    double newtonRaphson(double y, double x0, const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        const int maxIter = 100;
        const double tol = 1e-6;
        double xn = x0;

        for (int i = 0; i < maxIter; ++i) {
            double fx = evaluatePolynomial(xn, x, f, isForward) - y;
            double dfx = evaluateDerivative(xn, x, f, isForward);
            if (std::abs(dfx) < 1e-10) {
                return std::nan(""); // Avoid division by zero
            }
            double xn1 = xn - fx / dfx;
            if (std::abs(xn1 - xn) < tol) {
                return xn1;
            }
            xn = xn1;
        }
        return std::nan(""); // No convergence
    }

private slots:
    void createTable() {
        bool ok;
        int n = QInputDialog::getInt(this, "Number of Points", "Enter number of points:", 2, 2, 20, 1, &ok);
        if (!ok) return;

        tableWidget->clear();
        tableWidget->setRowCount(n);
        tableWidget->setColumnCount(n + 1);

        QStringList headers;
        headers << "X" << "f(X)";
        for (int i = 1; i <= n - 1; ++i) {
            headers << QString::number(i) + " Diff";
        }
        tableWidget->setHorizontalHeaderLabels(headers);
        equationLabel->setText("Polynomial Equation: ");
    }

    void computeDividedDifferences() {
        int n = tableWidget->rowCount();
        if (n == 0) return;

        QVector<QVector<double>> f(n, QVector<double>(n, 0));
        QVector<double> x(n);

        for (int i = 0; i < n; ++i) {
            x[i] = evaluateExpression(tableWidget->item(i, 0) ? tableWidget->item(i, 0)->text() : "0");
            f[0][i] = evaluateExpression(tableWidget->item(i, 1) ? tableWidget->item(i, 1)->text() : "0");
        }

        bool isForward = (methodCombo->currentText() == "Forward Interpolation");
        if (!isForward) {
            // Reverse x and f[0] for backward interpolation
            std::reverse(x.begin(), x.end());
            std::reverse(f[0].begin(), f[0].end());
        }

        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                f[i][j] = (f[i-1][j+1] - f[i-1][j]) / (x[j+i] - x[j]);
                tableWidget->setItem(j, i+1, new QTableWidgetItem(QString::number(f[i][j])));
            }
        }

        QString poly = generatePolynomial(x, f, isForward);
        equationLabel->setText("Polynomial Equation: " + poly);
    }

    void calculateValueAtX() {
        int n = tableWidget->rowCount();
        if (n == 0) return;

        QVector<QVector<double>> f(n, QVector<double>(n, 0));
        QVector<double> x(n);

        for (int i = 0; i < n; ++i) {
            x[i] = evaluateExpression(tableWidget->item(i, 0) ? tableWidget->item(i, 0)->text() : "0");
            f[0][i] = evaluateExpression(tableWidget->item(i, 1) ? tableWidget->item(i, 1)->text() : "0");
        }

        bool isForward = (methodCombo->currentText() == "Forward Interpolation");
        if (!isForward) {
            std::reverse(x.begin(), x.end());
            std::reverse(f[0].begin(), f[0].end());
        }

        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                f[i][j] = (f[i-1][j+1] - f[i-1][j]) / (x[j+i] - x[j]);
            }
        }

        bool ok;
        double XX = QInputDialog::getDouble(this, "Value of X", "Enter the value of X:", 0, -1e9, 1e9, 6, &ok);
        if (!ok) return;

        double P = evaluatePolynomial(XX, x, f, isForward);
        QMessageBox::information(this, "Result", "Interpolated value at x = " + QString::number(XX) + " is: " + QString::number(P));
    }

    void calculateXFromY() {
        int n = tableWidget->rowCount();
        if (n == 0) return;

        QVector<QVector<double>> f(n, QVector<double>(n, 0));
        QVector<double> x(n);

        for (int i = 0; i < n; ++i) {
            x[i] = evaluateExpression(tableWidget->item(i, 0) ? tableWidget->item(i, 0)->text() : "0");
            f[0][i] = evaluateExpression(tableWidget->item(i, 1) ? tableWidget->item(i, 1)->text() : "0");
        }

        bool isForward = (methodCombo->currentText() == "Forward Interpolation");
        if (!isForward) {
            std::reverse(x.begin(), x.end());
            std::reverse(f[0].begin(), f[0].end());
        }

        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                f[i][j] = (f[i-1][j+1] - f[i-1][j]) / (x[j+i] - x[j]);
            }
        }

        bool ok;
        double y = QInputDialog::getDouble(this, "Value of Y", "Enter the value of Y:", 0, -1e9, 1e9, 6, &ok);
        if (!ok) return;

        double x0 = QInputDialog::getDouble(this, "Initial Guess for X", "Enter initial guess for X:", 0, -1e9, 1e9, 6, &ok);
        if (!ok) return;

        double result = newtonRaphson(y, x0, x, f, isForward);
        if (std::isnan(result)) {
            QMessageBox::warning(this, "Error", "Newton-Raphson did not converge. Try a different initial guess.");
        } else {
            QMessageBox::information(this, "Result", "Value of x for y = " + QString::number(y) + " is: " + QString::number(result));
        }
    }

    void clearTable() {
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this, "Clear Table", "Are you sure you want to clear the table?",
                                      QMessageBox::Yes | QMessageBox::No);
        if (reply == QMessageBox::Yes) {
            tableWidget->clearContents();
            equationLabel->setText("Polynomial Equation: ");
        }
    }
};

#include "main.moc"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    DividedDifferencesApp w;
    w.setWindowTitle("Newton's Methods");
    w.resize(800, 600);
    w.show();

    return a.exec();
}


/*
#include <QApplication>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QTableWidget>
#include <QHeaderView>
#include <QInputDialog>
#include <QMessageBox>
#include <QPainter>
#include <QLabel>
#include <QComboBox>
#include <muParser.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <vector>

class DividedDifferencesApp : public QWidget {
    Q_OBJECT

public:
    DividedDifferencesApp(QWidget *parent = nullptr) : QWidget(parent) {
        setupUI();
    }

protected:
    void paintEvent(QPaintEvent *event) override {
        QPainter painter(this);
        QLinearGradient gradient(0, 0, width(), height());
        gradient.setColorAt(0, QColor(255, 255, 255));
        gradient.setColorAt(1, QColor(200, 230, 255));
        painter.fillRect(rect(), gradient);
        QWidget::paintEvent(event);
    }

private:
    QTableWidget *tableWidget;
    QPushButton *createButton, *computeButton, *calcButton, *calcXButton, *clearButton;
    QComboBox *methodCombo;
    QLabel *equationLabel;

    void setupUI() {
        QVBoxLayout *layout = new QVBoxLayout(this);

        tableWidget = new QTableWidget(this);
        tableWidget->setEditTriggers(QAbstractItemView::AllEditTriggers);
        tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
        tableWidget->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);
        tableWidget->setStyleSheet("QTableWidget { background-color: #ffffff; }");

        createButton = new QPushButton("Create Table", this);
        computeButton = new QPushButton("Compute Divided Differences", this);
        calcButton = new QPushButton("Calculate Value at X", this);
        calcXButton = new QPushButton("Calculate X from Y", this);
        clearButton = new QPushButton("Clear Table", this);

        QString buttonStyle = "QPushButton { background-color: #6ca0dc; color: white; font-weight: bold; padding: 8px; border-radius: 8px; }"
                              "QPushButton:hover { background-color: #5b8bd0; }";
        createButton->setStyleSheet(buttonStyle);
        computeButton->setStyleSheet(buttonStyle);
        calcButton->setStyleSheet(buttonStyle);
        calcXButton->setStyleSheet(buttonStyle);
        clearButton->setStyleSheet(buttonStyle);

        methodCombo = new QComboBox(this);
        methodCombo->addItem("Forward Interpolation");
        methodCombo->addItem("Backward Interpolation");
        methodCombo->setStyleSheet("QComboBox { padding: 5px; }");

        equationLabel = new QLabel("Polynomial Equation: ", this);
        equationLabel->setStyleSheet("QLabel { font-size: 14px; padding: 10px; }");

        QHBoxLayout *methodLayout = new QHBoxLayout();
        methodLayout->addWidget(new QLabel("Interpolation Method:", this));
        methodLayout->addWidget(methodCombo);

        layout->addWidget(tableWidget);
        layout->addLayout(methodLayout);
        layout->addWidget(equationLabel);
        layout->addWidget(createButton);
        layout->addWidget(computeButton);
        layout->addWidget(calcButton);
        layout->addWidget(calcXButton);
        layout->addWidget(clearButton);

        connect(createButton, &QPushButton::clicked, this, &DividedDifferencesApp::createTable);
        connect(computeButton, &QPushButton::clicked, this, &DividedDifferencesApp::computeDividedDifferences);
        connect(calcButton, &QPushButton::clicked, this, &DividedDifferencesApp::calculateValueAtX);
        connect(calcXButton, &QPushButton::clicked, this, &DividedDifferencesApp::calculateXFromY);
        connect(clearButton, &QPushButton::clicked, this, &DividedDifferencesApp::clearTable);
    }

    double evaluateExpression(const QString &expression) {
        try {
            mu::Parser parser;
            parser.DefineConst(L"pi", M_PI);
            parser.DefineConst(L"e", M_E);
            parser.SetExpr(expression.toStdWString());
            return parser.Eval();
        } catch (mu::Parser::exception_type &e) {
            QMessageBox::warning(this, "Error", QString::fromStdWString(e.GetMsg()));
            return 0.0;
        }
    }

    // Helper function to expand Newton polynomial to standard form
    std::vector<double> expandPolynomial(const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        int n = x.size();
        std::vector<double> coeffs(n, 0.0); // Coefficients of x^0, x^1, ..., x^(n-1)

        for (int i = 0; i < n; ++i) {
            double coef = f[i][0];
            if (std::abs(coef) < 1e-10) continue;

            // Temporary polynomial for (x-x_0)(x-x_1)...(x-x_{i-1})
            std::vector<double> temp_poly = {1.0}; // Start with 1
            for (int j = 0; j < i; ++j) {
                int idx = isForward ? j : (n - 1 - j);
                // Multiply temp_poly by (x - x[idx])
                std::vector<double> new_poly(temp_poly.size() + 1, 0.0);
                for (size_t k = 0; k < temp_poly.size(); ++k) {
                    new_poly[k] += temp_poly[k]; // x * temp_poly
                    new_poly[k + 1] += -x[idx] * temp_poly[k]; // -x[idx] * temp_poly
                }
                temp_poly = new_poly;
            }

            // Add coef * temp_poly to coeffs
            for (size_t k = 0; k < temp_poly.size(); ++k) {
                coeffs[k] += coef * temp_poly[k];
            }
        }

        return coeffs;
    }

    QString generatePolynomial(const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        std::vector<double> coeffs = expandPolynomial(x, f, isForward);
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);

        ss << "P(x) = ";
        bool firstTerm = true;
        for (int i = coeffs.size() - 1; i >= 0; --i) {
            double coef = coeffs[i];
            if (std::abs(coef) < 1e-10) continue;

            if (!firstTerm) {
                ss << (coef > 0 ? " + " : " - ");
            } else if (coef < 0) {
                ss << "-";
            }
            ss << std::abs(coef);
            if (i > 0) {
                ss << "x";
                if (i > 1) ss << "^" << i;
            }
            firstTerm = false;
        }
        if (firstTerm) ss << "0"; // Zero polynomial

        return QString::fromStdString(ss.str());
    }

    double evaluatePolynomial(double XX, const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        int n = x.size();
        double P = 0;
        for (int i = 0; i < n; ++i) {
            double k = 1;
            for (int j = 0; j < i; ++j) {
                int idx = isForward ? j : (n - 1 - j);
                k *= (XX - x[idx]);
            }
            P += k * f[i][0];
        }
        return P;
    }

    double evaluateDerivative(double XX, const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        int n = x.size();
        double deriv = 0;
        for (int i = 1; i < n; ++i) {
            double term = f[i][0];
            for (int j = 0; j < i; ++j) {
                double prod = 1;
                for (int k = 0; k < i; ++k) {
                    if (k == j) continue;
                    int idx = isForward ? k : (n - 1 - k);
                    prod *= (XX - x[idx]);
                }
                deriv += term * prod;
            }
        }
        return deriv;
    }

    double newtonRaphson(double y, double x0, const QVector<double> &x, const QVector<QVector<double>> &f, bool isForward) {
        const int maxIter = 100;
        const double tol = 1e-6;
        double xn = x0;

        for (int i = 0; i < maxIter; ++i) {
            double fx = evaluatePolynomial(xn, x, f, isForward) - y;
            double dfx = evaluateDerivative(xn, x, f, isForward);
            if (std::abs(dfx) < 1e-10) {
                return std::nan(""); // Avoid division by zero
            }
            double xn1 = xn - fx / dfx;
            if (std::abs(xn1 - xn) < tol) {
                return xn1;
            }
            xn = xn1;
        }
        return std::nan(""); // No convergence
    }

private slots:
    void createTable() {
        bool ok;
        int n = QInputDialog::getInt(this, "Number of Points", "Enter number of points:", 2, 2, 20, 1, &ok);
        if (!ok) return;

        tableWidget->clear();
        tableWidget->setRowCount(n);
        tableWidget->setColumnCount(n + 1);

        QStringList headers;
        headers << "X" << "f(X)";
        for (int i = 1; i <= n - 1; ++i) {
            headers << QString::number(i) + " Diff";
        }
        tableWidget->setHorizontalHeaderLabels(headers);
        equationLabel->setText("Polynomial Equation: ");
    }

    void computeDividedDifferences() {
        int n = tableWidget->rowCount();
        if (n == 0) return;

        QVector<QVector<double>> f(n, QVector<double>(n, 0));
        QVector<double> x(n);

        for (int i = 0; i < n; ++i) {
            x[i] = evaluateExpression(tableWidget->item(i, 0) ? tableWidget->item(i, 0)->text() : "0");
            f[0][i] = evaluateExpression(tableWidget->item(i, 1) ? tableWidget->item(i, 1)->text() : "0");
        }

        bool isForward = (methodCombo->currentText() == "Forward Interpolation");
        if (!isForward) {
            std::reverse(x.begin(), x.end());
            std::reverse(f[0].begin(), f[0].end());
        }

        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                f[i][j] = (f[i-1][j+1] - f[i-1][j]) / (x[j+i] - x[j]);
                tableWidget->setItem(j, i+1, new QTableWidgetItem(QString::number(f[i][j])));
            }
        }

        QString poly = generatePolynomial(x, f, isForward);
        equationLabel->setText("Polynomial Equation: " + poly);
    }

    void calculateValueAtX() {
        int n = tableWidget->rowCount();
        if (n == 0) return;

        QVector<QVector<double>> f(n, QVector<double>(n, 0));
        QVector<double> x(n);

        for (int i = 0; i < n; ++i) {
            x[i] = evaluateExpression(tableWidget->item(i, 0) ? tableWidget->item(i, 0)->text() : "0");
            f[0][i] = evaluateExpression(tableWidget->item(i, 1) ? tableWidget->item(i, 1)->text() : "0");
        }

        bool isForward = (methodCombo->currentText() == "Forward Interpolation");
        if (!isForward) {
            std::reverse(x.begin(), x.end());
            std::reverse(f[0].begin(), f[0].end());
        }

        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                f[i][j] = (f[i-1][j+1] - f[i-1][j]) / (x[j+i] - x[j]);
            }
        }

        bool ok;
        double XX = QInputDialog::getDouble(this, "Value of X", "Enter the value of X:", 0, -1e9, 1e9, 6, &ok);
        if (!ok) return;

        double P = evaluatePolynomial(XX, x, f, isForward);
        QMessageBox::information(this, "Result", "Interpolated value at x = " + QString::number(XX) + " is: " + QString::number(P));
    }

    void calculateXFromY() {
        int n = tableWidget->rowCount();
        if (n == 0) return;

        QVector<QVector<double>> f(n, QVector<double>(n, 0));
        QVector<double> x(n);

        for (int i = 0; i < n; ++i) {
            x[i] = evaluateExpression(tableWidget->item(i, 0) ? tableWidget->item(i, 0)->text() : "0");
            f[0][i] = evaluateExpression(tableWidget->item(i, 1) ? tableWidget->item(i, 1)->text() : "0");
        }

        bool isForward = (methodCombo->currentText() == "Forward Interpolation");
        if (!isForward) {
            std::reverse(x.begin(), x.end());
            std::reverse(f[0].begin(), f[0].end());
        }

        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                f[i][j] = (f[i-1][j+1] - f[i-1][j]) / (x[j+i] - x[j]);
            }
        }

        bool ok;
        double y = QInputDialog::getDouble(this, "Value of Y", "Enter the value of Y:", 0, -1e9, 1e9, 6, &ok);
        if (!ok) return;

        double x0 = QInputDialog::getDouble(this, "Initial Guess for X", "Enter initial guess for X:", 0, -1e9, 1e9, 6, &ok);
        if (!ok) return;

        double result = newtonRaphson(y, x0, x, f, isForward);
        if (std::isnan(result)) {
            QMessageBox::warning(this, "Error", "Newton-Raphson did not converge. Try a different initial guess.");
        } else {
            QMessageBox::information(this, "Result", "Value of x for y = " + QString::number(y) + " is: " + QString::number(result));
        }
    }

    void clearTable() {
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this, "Clear Table", "Are you sure you want to clear the table?",
                                      QMessageBox::Yes | QMessageBox::No);
        if (reply == QMessageBox::Yes) {
            tableWidget->clearContents();
            equationLabel->setText("Polynomial Equation: ");
        }
    }
};

#include "main.moc"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    DividedDifferencesApp w;
    w.setWindowTitle("Newton's Methods");
    w.resize(800, 600);
    w.show();

    return a.exec();
}
*/





