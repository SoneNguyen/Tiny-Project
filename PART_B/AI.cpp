#include "../PART_A/Vector.h"
#include "../PART_A/Matrix.h"
#include "../PART_A/LinearSystem.h"
#include "../PART_A/Tikhonov.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>
#include <iomanip>

struct CPUData {
    std::string vendor;
    std::string model;
    double myct;    // machine cycle time in nanoseconds
    double mmin;    // minimum main memory in kilobytes
    double mmax;    // maximum main memory in kilobytes
    double cach;    // cache memory in kilobytes
    double chmin;   // minimum channels in units
    double chmax;   // maximum channels in units
    double prp;     // published relative performance (target)
    double erp;     // estimated relative performance
};

    class CPUPerformanceRegression {
private:
    std::vector<CPUData> data;
    std::vector<CPUData> trainData;
    std::vector<CPUData> testData;
    Vector parameters; // model parameters

    // Helper: safe std deviation (avoid zero)
    double safeStd(double stdVal) {
        return stdVal < 1e-12 ? 1.0 : stdVal;
    }

public:
    // Load data from CSV file
    bool loadData(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string item;
            CPUData cpu;

            try {
                std::getline(ss, cpu.vendor, ',');
                std::getline(ss, cpu.model, ',');

                std::getline(ss, item, ',');
                cpu.myct = std::stod(item);

                std::getline(ss, item, ',');
                cpu.mmin = std::stod(item);

                std::getline(ss, item, ',');
                cpu.mmax = std::stod(item);

                std::getline(ss, item, ',');
                cpu.cach = std::stod(item);

                std::getline(ss, item, ',');
                cpu.chmin = std::stod(item);

                std::getline(ss, item, ',');
                cpu.chmax = std::stod(item);

                std::getline(ss, item, ',');
                cpu.prp = std::stod(item);

                std::getline(ss, item, ',');
                cpu.erp = std::stod(item);

                data.push_back(cpu);
            } catch (const std::exception& e) {
                std::cerr << "Error parsing line: " << line << std::endl;
                continue;
            }
        }

        file.close();
        std::cout << "Loaded " << data.size() << " data points" << std::endl;
        return true;
    }

    //Feature Engineering
    Vector getFeatures(const CPUData& cpu) {
        Vector features(10);

        // 6 base value
        features(1) = cpu.myct;
        features(2) = cpu.mmin;
        features(3) = cpu.mmax;
        features(4) = cpu.cach;
        features(5) = cpu.chmin;
        features(6) = cpu.chmax;

        // 4 new var
        features(7) = (cpu.mmax != 0) ? (cpu.mmin / cpu.mmax) : 0;  // ratio mmin/mmax
        features(8) = cpu.myct * cpu.cach;                         // myct*cach
        features(9) = cpu.myct * cpu.myct;                         // myct^2
        features(10) = cpu.cach * cpu.cach;                        // cach^2

        return features;
    }
    
    // Split data into training (80%) and testing (20%)
    void splitData(double trainRatio = 0.8) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(data.begin(), data.end(), g);

        size_t trainSize = static_cast<size_t>(data.size() * trainRatio);

        trainData.assign(data.begin(), data.begin() + trainSize);
        testData.assign(data.begin() + trainSize, data.end());

        std::cout << "Training set: " << trainData.size() << " samples" << std::endl;
        std::cout << "Testing set: " << testData.size() << " samples" << std::endl;
    }

    // Normalize features using training set stats
    void     normalizeFeatures() {
        if (trainData.empty()) return;

        // Calculate mean
        double meanMyct = 0, meanMmin = 0, meanMmax = 0, meanCach = 0, meanChmin = 0, meanChmax = 0;

        for (const auto& cpu : trainData) {
            meanMyct += cpu.myct;
            meanMmin += cpu.mmin;
            meanMmax += cpu.mmax;
            meanCach += cpu.cach;
            meanChmin += cpu.chmin;
            meanChmax += cpu.chmax;
        }

        size_t n = trainData.size();
        meanMyct /= n; meanMmin /= n; meanMmax /= n;
        meanCach /= n; meanChmin /= n; meanChmax /= n;

        // Calculate std dev
        double stdMyct = 0, stdMmin = 0, stdMmax = 0, stdCach = 0, stdChmin = 0, stdChmax = 0;

        for (const auto& cpu : trainData) {
            stdMyct += (cpu.myct - meanMyct) * (cpu.myct - meanMyct);
            stdMmin += (cpu.mmin - meanMmin) * (cpu.mmin - meanMmin);
            stdMmax += (cpu.mmax - meanMmax) * (cpu.mmax - meanMmax);
            stdCach += (cpu.cach - meanCach) * (cpu.cach - meanCach);
            stdChmin += (cpu.chmin - meanChmin) * (cpu.chmin - meanChmin);
            stdChmax += (cpu.chmax - meanChmax) * (cpu.chmax - meanChmax);
        }

        stdMyct = safeStd(std::sqrt(stdMyct / (n - 1)));
        stdMmin = safeStd(std::sqrt(stdMmin / (n - 1)));
        stdMmax = safeStd(std::sqrt(stdMmax / (n - 1)));
        stdCach = safeStd(std::sqrt(stdCach / (n - 1)));
        stdChmin = safeStd(std::sqrt(stdChmin / (n - 1)));
        stdChmax = safeStd(std::sqrt(stdChmax / (n - 1)));

        // Normalize training data
        for (auto& cpu : trainData) {
            cpu.myct = (cpu.myct - meanMyct) / stdMyct;
            cpu.mmin = (cpu.mmin - meanMmin) / stdMmin;
            cpu.mmax = (cpu.mmax - meanMmax) / stdMmax;
            cpu.cach = (cpu.cach - meanCach) / stdCach;
            cpu.chmin = (cpu.chmin - meanChmin) / stdChmin;
            cpu.chmax = (cpu.chmax - meanChmax) / stdChmax;
        }

        // Normalize testing data using train stats
        for (auto& cpu : testData) {
            cpu.myct = (cpu.myct - meanMyct) / stdMyct;
            cpu.mmin = (cpu.mmin - meanMmin) / stdMmin;
            cpu.mmax = (cpu.mmax - meanMmax) / stdMmax;
            cpu.cach = (cpu.cach - meanCach) / stdCach;
            cpu.chmin = (cpu.chmin - meanChmin) / stdChmin;
            cpu.chmax = (cpu.chmax - meanChmax) / stdChmax;
        }

        std::cout << "Features normalized" << std::endl;
    }

    // Calculate RMSE between predicted and actual
    double calculateRMSE(const Vector& predictions, const Vector& actual) {
        double sum = 0.0;
        int n = predictions.Size();

        for (int i = 1; i <= n; i++) {
            double diff = predictions(i) - actual(i);
            sum += diff * diff;
        }

        return std::sqrt(sum / n);
    }

    // Grid search for best alpha (regularization parameter)
    double gridSearchAlpha(const std::vector<double>& alphas) {
        if (trainData.empty()) {
            std::cerr << "No training data for Grid Search" << std::endl;
            return 0.01;
        }

        double bestAlpha = alphas[0];
        double bestRmse = 1e9;
        Vector bestParams;

        int numSamples = static_cast<int>(trainData.size());
        int numFeatures = 10;  //10 features

        Vector actual(numSamples);
        for (int i = 0; i < numSamples; i++) {
            actual(i + 1) = trainData[i].prp;
        }

        for (double alpha : alphas) {
            Matrix A(numSamples, numFeatures + 1); // +1 bias
            Vector b(numSamples);

            for (int i = 0; i < numSamples; i++) {
                Vector feats = getFeatures(trainData[i]);
                for (int j = 1; j <= numFeatures; j++) {
                    A(i + 1, j) = feats(j);
                }
                A(i + 1, numFeatures + 1) = 1.0; // bias
                b(i + 1) = trainData[i].prp;
            }

            try {
                TikhonovSolver tikSolver(alpha);
                Vector params = tikSolver.Solve(A, b);

                // Dự đoán trên train
                Vector predictions(numSamples);
                for (int i = 0; i < numSamples; i++) {
                    Vector feats = getFeatures(trainData[i]);
                    double pred = 0;
                    for (int j = 1; j <= numFeatures; j++) {
                        pred += params(j) * feats(j);
                    }
                    pred += params(numFeatures + 1);
                    predictions(i + 1) = pred;
                }

                double rmse = calculateRMSE(predictions, actual);
                std::cout << "Alpha = " << alpha << ", RMSE = " << rmse << std::endl;

                if (rmse < bestRmse) {
                    bestRmse = rmse;
                    bestAlpha = alpha;
                    bestParams = params;
                }
            } catch (const std::exception& e) {
                std::cerr << "Error with alpha = " << alpha << ": " << e.what() << std::endl;
            }
        }

        parameters = bestParams;
        std::cout << "Best alpha found: " << bestAlpha << " with RMSE = " << bestRmse << std::endl;
        return bestAlpha;
    }

    // Train model
    void train() {
        if (trainData.empty()) {
            std::cerr << "No training data available" << std::endl;
            return;
        }

        std::vector<double> alphas = {0.0001, 0.001, 0.01, 0.1, 1.0, 10.0};
        double bestAlpha = gridSearchAlpha(alphas);

        int numSamples = static_cast<int>(trainData.size());
        int numFeatures = 10;

        Matrix A(numSamples, numFeatures + 1);
        Vector b(numSamples);

        for (int i = 0; i < numSamples; i++) {
            Vector feats = getFeatures(trainData[i]);
            for (int j = 1; j <= numFeatures; j++) {
                A(i + 1, j) = feats(j);
            }
            A(i + 1, numFeatures + 1) = 1.0; // bias
            b(i + 1) = trainData[i].prp;
        }

        try {
            TikhonovSolver tikSolver(bestAlpha);
            parameters = tikSolver.Solve(A, b);

            std::cout << "Training completed with alpha = " << bestAlpha << std::endl;
            std::cout << "Model parameters:" << std::endl;
            for (int i = 1; i <= numFeatures; i++) {
                std::cout << "x" << i << ": " << parameters(i) << std::endl;
            }
            std::cout << "Bias (x" << (numFeatures + 1) << "): " << parameters(numFeatures + 1) << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error during training: " << e.what() << std::endl;
        }
    }
    // Predict given dataset
    Vector predict(const std::vector<CPUData>& dataset) {
        int numFeatures = 10;
        Vector predictions(static_cast<int>(dataset.size()));

        for (size_t i = 0; i < dataset.size(); i++) {
            Vector feats = getFeatures(dataset[i]);
            double prediction = 0;
            for (int j = 1; j <= numFeatures; j++) {
                prediction += parameters(j) * feats(j);
            }
            prediction += parameters(numFeatures + 1);
            prediction = std::max(0.0, prediction);
            predictions(i + 1) = prediction;
        }

        return predictions;
    }

    // Evaluate model performance
    void evaluate() {
        if (testData.empty()) {
            std::cerr << "No test data available" << std::endl;
            return;
        }
        
        // Get predictions for test set
        Vector predictions = predict(testData);
        
        // Create actual values vector
        Vector actual(testData.size());
        for (size_t i = 0; i < testData.size(); i++) {
            actual(i+1) = testData[i].prp;
        }
        
        // Calculate RMSE
        double rmse = calculateRMSE(predictions, actual);
        
        std::cout << "\n=== Model Evaluation ===" << std::endl;
        std::cout << "Test set RMSE: " << rmse << std::endl;
        
        // Show some sample predictions vs actual
        std::cout << "\nSample predictions vs actual:" << std::endl;
        std::cout << "Predicted\tActual\t\tDifference" << std::endl;
        int sampleSize = std::min(10, static_cast<int>(testData.size()));
        for (int i = 1; i <= sampleSize; i++) {
            double pred = predictions(i);
            double act = actual(i);
            std::cout << std::fixed << std::setprecision(2);
            std::cout << pred << "\t\t" << act << "\t\t" << (pred - act) << std::endl;
        }
    }
    
    // Run the complete pipeline
    void run(const std::string& filename) {
        std::cout << "=== CPU Performance Linear Regression ===" << std::endl;
        
        // Load data
        if (!loadData(filename)) {
            return;
        }
        
        // Split data
        splitData(0.8);
        
        // Normalize features for better numerical stability
        normalizeFeatures();
        
        // Train model
        train();
        
        // Evaluate model
        evaluate();
    }
};

int main(void) {
    CPUPerformanceRegression regression;
    
    // Replace with actual path to your dataset
    std::string dataFile = "machine.data"; // UCI dataset filename
    
    regression.run(dataFile);
    
    return 0;
}