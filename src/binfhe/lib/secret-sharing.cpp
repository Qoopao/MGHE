#include "secret-sharing.h"
#include <bits/stdint-uintn.h>
#include <vector>
#include <filesystem>
#include <string>
#include "math/hal/nativeintbackend.h"
#include "utils/inttypes.h"




namespace lbcrypto {

// 生成 triples
void AdditiveSecretSharing::TriplesGen(BinFHEContext& cc, int groupidx, int numofparties, int numofTriples) {
    if (numofparties <= 0) {
        throw std::invalid_argument("numofparties must be positive");
    }

    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    std::string group_dir = user_dir + "/group_" + std::to_string(groupidx);

    auto check = [](const auto& path, const auto& msg) {
        if (!std::filesystem::exists(path)) throw std::runtime_error(msg + path);
    };

    // 检查必要文件
    check(user_dir + "/rlweCRP", "rlweCRP file not found: ");
    check(group_dir + "/jointrlwepk", "jointrlwepk file not found: ");

    // 检查所有party文件
    for (int i = 0; i < numofparties; ++i) {
        std::string prefix = group_dir + "/party_" + std::to_string(i);
        check(prefix + "/rlwesk", "rlwesk file not found: ");
        check(prefix + "/rlwepk", "rlwepk file not found: ");
    }


    std::vector<std::vector<std::vector<NativePoly>>> triples;
    triples.resize(numofTriples);

    for(int i=0;i<numofTriples;++i){

        std::vector<std::vector<NativePoly>> triple(3,std::vector<NativePoly>(numofparties));
        auto params = cc.GetParams()->GetVectorNTRUParams();
        TernaryUniformGeneratorImpl<NativeVector> tug;
        std::vector<RLWECiphertext> C(numofparties);
        for(int i=0;i<numofparties;++i){       //参与方P_i生成随机多项式a_i，加密
            NativePoly a_poly(tug,params->GetPolyParams(),Format::COEFFICIENT);
            a_poly.SetFormat(Format::EVALUATION);
            triple[0][i] = a_poly;
            C[i] = cc.AHE_Encrypt(cc.GetParams(),a_poly,0);
        }

        //P_0计算密文和，得到a的密文
        RLWECiphertext a_sum = C[0];
        for(int i=1;i<numofparties;++i){
            a_sum = cc.AHE_Add(a_sum,C[i]);
        }


        //各方计算E_i=Enc(0)，随机选取b_i，计算C_i'=Add(ScalarMult(E_i,b_i),a_sum)
        std::vector<RLWECiphertext> C_prime(numofparties);
        for(int i=0;i<numofparties;++i){      
            NativePoly b_poly(tug,params->GetPolyParams(),Format::COEFFICIENT);
            b_poly.SetFormat(Format::EVALUATION);
            triple[1][i] = b_poly;
            NativePoly zero = NativePoly(params->GetPolyParams(),Format::EVALUATION,true);
            RLWECiphertext E_i = cc.AHE_Encrypt(cc.GetParams(),zero,0);
            C_prime[i] = cc.AHE_Add(cc.AHE_ScalarMult(a_sum,b_poly),E_i);
        }

        //除P_0，其他方随机选取c_i，计算C_i''=Enc(-c_i)
        std::vector<RLWECiphertext> C_prime_prime(numofparties); 
        for(int i=1;i<numofparties;++i){
            NativePoly c_poly(tug,params->GetPolyParams(),Format::COEFFICIENT);
            c_poly.SetFormat(Format::EVALUATION);
            triple[2][i] = c_poly;
            auto c_poly_negate = c_poly.Negate();
            C_prime_prime[i] = cc.AHE_Encrypt(cc.GetParams(),c_poly_negate,0);
        }
        //P_0计算所有C_i'的和再加上所有C_i''的和，得到ab-\sum c_i的密文，然后解密得到c_0
        RLWECiphertext ab_minus_sum_c = C_prime[0];
        for(int i=1;i<numofparties;++i){
            ab_minus_sum_c = cc.AHE_Add(ab_minus_sum_c,C_prime[i]);
            ab_minus_sum_c = cc.AHE_Add(ab_minus_sum_c,C_prime_prime[i]);
        }

        auto c_0 = cc.AHE_Decrypt(cc.GetParams(),ab_minus_sum_c,0,numofparties);
        c_0.SetFormat(Format::EVALUATION);
        triple[2][0] = c_0;

        
        triples[i] = triple;
        
    }
    
    this->setTriples(triples);

}


std::vector<NativePoly> AdditiveSecretSharing::Split(NativePoly& s, int k) {
    if (k <= 0) {
        throw std::invalid_argument("k must be positive");
    }
    if(s.GetFormat() != Format::EVALUATION) {
        throw std::invalid_argument("s must be in EVALUATION format");
    }
    
    std::vector<NativePoly> shares(k, NativePoly(vectorNTRUParams->GetPolyParams(), Format::EVALUATION, true));
    NativePoly sum(vectorNTRUParams->GetPolyParams(), Format::EVALUATION, true); // 初始化为0多项式
    
    // 生成k-1个随机分片
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(vectorNTRUParams->GetQ());
    for (int i = 0; i < k - 1; ++i) {
        NativePoly randomPoly(dug, vectorNTRUParams->GetPolyParams(), Format::EVALUATION);
        shares[i] =  randomPoly;
        sum = sum + shares[i];
    }
    
    // 计算最后一个分片：s - sum
    shares[k - 1] = s - sum;
    
    return shares;
}

NativePoly AdditiveSecretSharing::Recover(std::vector<NativePoly>& shares) {
    if (shares.empty()) {
        throw std::invalid_argument("shares cannot be empty");
    }
    NativePoly result(vectorNTRUParams->GetPolyParams(), Format::EVALUATION, true); // 初始化为0多项式
    for (auto& share : shares) {
        result = result + share;
    }
    
    return result;
}

std::vector<NativePoly> AdditiveSecretSharing::Add(
    std::vector<NativePoly>& x_shares,    
    std::vector<NativePoly>& y_shares) {    
    
    if (x_shares.size() != y_shares.size()) {
        throw std::invalid_argument("x_shares and y_shares must have same size");
    }
    
    std::vector<NativePoly> z_shares(x_shares.size(), NativePoly(vectorNTRUParams->GetPolyParams(), Format::EVALUATION, true));
    z_shares.reserve(x_shares.size());
    
    for (size_t i = 0; i < x_shares.size(); ++i) {
        z_shares[i] = x_shares[i] + y_shares[i];
    }
    
    return z_shares;
}

std::vector<NativePoly> AdditiveSecretSharing::ScalarMult(
    NativePoly& alpha,
    std::vector<NativePoly>& x_shares) {
    
    if(alpha.GetFormat() != Format::EVALUATION) {
        throw std::invalid_argument("alpha must be in EVALUATION format");
    }
        
    std::vector<NativePoly> z_shares(x_shares.size(), NativePoly(vectorNTRUParams->GetPolyParams(), Format::EVALUATION, true));
    z_shares.reserve(x_shares.size());
    
    for (size_t i = 0; i < x_shares.size(); ++i) {
        z_shares[i] = alpha * x_shares[i];
    }

    return z_shares;
}

std::vector<NativePoly> AdditiveSecretSharing::Mult(
    std::vector<NativePoly>& x_shares,
    std::vector<NativePoly>& y_shares,
    std::vector<NativePoly>& a_shares,
    std::vector<NativePoly>& b_shares,
    std::vector<NativePoly>& c_shares) {
    
    int k = x_shares.size();
    if (y_shares.size() != k || a_shares.size() != k || 
        b_shares.size() != k || c_shares.size() != k) {
        throw std::invalid_argument("All shares must have same size");
    }
    
    // 计算 [x]_i' = [x]_i - [a]_i 和 [y]_i' = [y]_i - [b]_i
    std::vector<NativePoly> x_prime_shares, y_prime_shares;
    for (int i = 0; i < k; ++i) {
        x_prime_shares.push_back(x_shares[i] - a_shares[i]);
        y_prime_shares.push_back(y_shares[i] - b_shares[i]);
    }
    
    // 恢复 x' 和 y'
    NativePoly x_prime = Recover(x_prime_shares);
    NativePoly y_prime = Recover(y_prime_shares);
    
    // 计算 [z]_i
    std::vector<NativePoly> z_shares(k, NativePoly(vectorNTRUParams->GetPolyParams(), Format::EVALUATION, true));
    
    // P1 的特殊计算
    z_shares[0] = x_prime * y_prime + c_shares[0] + 
                   y_prime * a_shares[0] + x_prime * b_shares[0];

    // 其他参与方的计算
    for (int i = 1; i < k; ++i) {
        z_shares[i] = c_shares[i] + y_prime * a_shares[i] + x_prime * b_shares[i];
    }

    return z_shares;
}

std::vector<NativePoly> AdditiveSecretSharing::Auto(
    std::vector<NativePoly>& x_shares,
    uint32_t& sigma) {
    
    std::vector<NativePoly> z_shares(x_shares.size(), NativePoly(vectorNTRUParams->GetPolyParams(), Format::COEFFICIENT, true));
    z_shares.reserve(x_shares.size());
    
    for (size_t i = 0; i < x_shares.size(); ++i) {
        z_shares[i] = x_shares[i].AutomorphismTransform(sigma);
    }
    
    return z_shares;
}

}