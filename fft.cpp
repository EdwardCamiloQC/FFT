#include <fft.h>
#include <stdio.h>
#include <iostream>

void calculateFFT(const double* signal, unsigned int length, std::complex<double>* transform){
    if(length>1){
        double par[int(length/2)], impar[int(length/2)];
        std::complex<double> oddTransform[int(length/2)], evenTransform[int(length/2)];
        oddValues(signal, length, par);
        evenValues(signal, length, impar);
        calculateFFT(par, sizeof(par)/sizeof(double), oddTransform);
        calculateFFT(impar, sizeof(impar)/sizeof(double), evenTransform);
        for(size_t k=0; k<length; k++){
            if(k<length/2){
                transform[k] = oddTransform[k] + evenTransform[k] * std::polar(1.0,-2.0*M_PI*k/length);
                transform[k+length/2] = oddTransform[k] - evenTransform[k] * std::polar(1.0,-2.0*M_PI*k/length);
            }
        }
    }if(length==1){
        *transform = *signal;
    }
}

void calculateFFT(const std::complex<double>* signal, unsigned int length, std::complex<double>* transform){
    if(length>1){
        std::complex<double> par[int(length/2)], impar[int(length/2)];
        std::complex<double> oddTransform[int(length/2)], evenTransform[int(length/2)];
        oddValues(signal, length, par);
        evenValues(signal, length, impar);
        calculateFFT(par, sizeof(par)/sizeof(std::complex<double>), oddTransform);
        calculateFFT(impar, sizeof(impar)/sizeof(std::complex<double>), evenTransform);
        for(size_t k=0; k<length; k++){
            if(k<length/2){
                transform[k] = oddTransform[k] + evenTransform[k] * std::polar(1.0,-2.0*M_PI*k/length);
                transform[k+length/2] = oddTransform[k] - evenTransform[k] * std::polar(1.0,-2.0*M_PI*k/length);
            }
        }
    }if(length==1){
        *transform = *signal;
    }
}

void calculateFFT(const float* signal, unsigned int length, std::complex<float>* transform){
    if(length>1){
        float par[int(length/2)], impar[int(length/2)];
        std::complex<float> oddTransform[int(length/2)], evenTransform[int(length/2)];
        oddValues(signal, length, par);
        evenValues(signal, length, impar);
        calculateFFT(par, sizeof(par)/sizeof(float), oddTransform);
        calculateFFT(impar, sizeof(impar)/sizeof(float), evenTransform);
        for(size_t k=0; k<length; k++){
            if(k<length/2){
                transform[k] = oddTransform[k] + evenTransform[k] * std::polar<float>(1.0,-2.0*M_PI*k/length);
                transform[k+length/2] = oddTransform[k] - evenTransform[k] * std::polar<float>(1.0,-2.0*M_PI*k/length);
            }
        }
    }if(length==1){
        *transform = *signal;
    }
}

void calculateFFT(const std::complex<float>* signal, unsigned int length, std::complex<float>* transform){
    if(length>1){
        std::complex<float> par[int(length/2)], impar[int(length/2)];
        std::complex<float> oddTransform[int(length/2)], evenTransform[int(length/2)];
        oddValues(signal, length, par);
        evenValues(signal, length, impar);
        calculateFFT(par, sizeof(par)/sizeof(std::complex<float>), oddTransform);
        calculateFFT(impar, sizeof(impar)/sizeof(std::complex<float>), evenTransform);
        for(size_t k=0; k<length; k++){
            if(k<length/2){
                transform[k] = oddTransform[k] + evenTransform[k] * std::polar<float>(1.0,-2.0*M_PI*k/length);
                transform[k+length/2] = oddTransform[k] - evenTransform[k] * std::polar<float>(1.0,-2.0*M_PI*k/length);
            }
        }
    }if(length==1){
        *transform = *signal;
    }
}

void calculateModule(const std::complex<double>* transform, unsigned int n, double* module){
    for(size_t k=0; k<n; k++){
        module[k] = std::abs(transform[k]);
    }
}

void seeTransform(const std::complex<double>* transform, unsigned int n){
    for(size_t k=0; k<n; k++){
        std::cout << transform[k];
    }
}

void seeTransform(const std::complex<float>* transform, unsigned int n){
    for(size_t k=0; k<n; k++){
        std::cout << transform[k];
    }
}
