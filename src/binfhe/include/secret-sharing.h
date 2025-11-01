#ifndef SECRET_SHARING_H
#define SECRET_SHARING_H

#include <cstdint>
#include <vector>
#include "binfhe-constants.h"
#include "binfhecontext.h"
#include "lattice/hal/lat-backend.h"
#include "math/hal/nativeintbackend.h"


namespace lbcrypto {

// 加性秘密共享类
class AdditiveSecretSharing {
private:
    std::shared_ptr<VectorNTRUCryptoParams> vectorNTRUParams;

    //从外到内分别是三元组的个数，3（a,b,c），参与方个数
    std::vector<std::vector<std::vector<NativePoly>>> Triples; 
    
    int k;

    int numofTriples;

public:
    // 构造函数
    AdditiveSecretSharing(const std::shared_ptr<VectorNTRUCryptoParams>& vectorNTRUParams, int k, int numofTriples = 1){
        this->vectorNTRUParams = vectorNTRUParams;
        this->k = k;
        this->numofTriples = numofTriples;
        Triples.resize(numofTriples);
        for (int i = 0; i < numofTriples; i++){
            Triples[i].resize(3);
            for (int j = 0; j < 3; j++){
                Triples[i][j].resize(k);
            }
        }
    }
    
    void setTriples(std::vector<std::vector<std::vector<NativePoly>>>& triples){
        Triples = triples;
    }
    
    std::vector<std::vector<std::vector<NativePoly>>> getTriples(){
        return Triples;
    }

    void setk(int k){
        this->k = k;
    }
    
    int getk(){
        return k;
    }

    void setnumofTriples(int numofTriples){
        this->numofTriples = numofTriples;
    }
    
    int getnumofTriples(){
        return numofTriples;
    }

    void TriplesGen(BinFHEContext& cc, int groupidx, int numofparties, int numofTriples = 1);

    // 分割秘密
    std::vector<NativePoly> Split(NativePoly& s, int k);
    
    // 恢复秘密
    NativePoly Recover(std::vector<NativePoly>& shares);
    
    // 加法
    std::vector<NativePoly> Add(
        std::vector<NativePoly>& x_shares,
        std::vector<NativePoly>& y_shares);
    
    // 标量乘法
    std::vector<NativePoly> ScalarMult(
        NativePoly& alpha,            
        std::vector<NativePoly>& x_shares);
    
    // 乘法（使用三元组）
    std::vector<NativePoly> Mult(
        std::vector<NativePoly>& x_shares,
        std::vector<NativePoly>& y_shares,
        std::vector<NativePoly>& a_shares,
        std::vector<NativePoly>& b_shares,
        std::vector<NativePoly>& c_shares);
    
    // 自同构
    std::vector<NativePoly> Auto(
        std::vector<NativePoly>& x_shares,
        uint32_t& sigma); 
};

#endif // SECRET_SHARING_H


} // namespace lbcrypto 