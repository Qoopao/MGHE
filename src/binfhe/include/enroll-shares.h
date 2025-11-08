#ifndef ENROLL_SHARES_H
#define ENROLL_SHARES_H

#include <cstdint>
#include <iostream>
#include <vector>
#include "binfhe-constants.h"
#include "binfhecontext.h"
#include "lattice/hal/lat-backend.h"
#include "math/hal/nativeintbackend.h"
#include "secret-sharing.h"
#include <filesystem>

namespace fs = std::filesystem;

namespace lbcrypto{

/**
 * @brief 存储参与方注册过程中的共享数据
 */
class EnrollShares {
private:
    // F的共享
    std::vector<lbcrypto::NativePoly> F_shares;
    // F逆的共享
    std::vector<lbcrypto::NativePoly> F_inv_shares;
    // 私钥逆的和的共享
    std::vector<std::vector<lbcrypto::NativePoly>> sk_sum_shares;
    // 私钥全部和的共享
    std::vector<lbcrypto::NativePoly> sk_inv_all_sum_shares;

public:
    EnrollShares() = default;
    
    /**
     * @brief 根据参与方数量初始化
     * @param numofparties 参与方数量
     * @param n 私钥长度
     */

    //根据论文内容，这里只能串行：先把所有参与方的私钥读进来，然后按论文的伪代码把生成过程放到同一个循环里面
    EnrollShares(BinFHEContext& cc, int groupid, int numofparties, uint32_t n) {
        F_shares.resize(numofparties);
        F_inv_shares.resize(numofparties);
        sk_sum_shares.resize(n);
        for(int partyidx=0;partyidx<numofparties;++partyidx){sk_sum_shares[partyidx].resize(numofparties);}
        sk_inv_all_sum_shares.resize(numofparties);
    
        auto params = cc.GetParams()->GetVectorNTRUParams();
        uint32_t N = params->GetN();
        NativeInteger Q = params->GetQ();
        AdditiveSecretSharing secretSharing(cc.GetParams()->GetVectorNTRUParams(), numofparties);
        secretSharing.TriplesGen(cc, groupid, numofparties);
        auto Triples = secretSharing.getTriples();
        auto a_shares = Triples[0][0];
        auto b_shares = Triples[0][1];
        auto c_shares = Triples[0][2];

        std::string groupDir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupid);

        //生成F与F^-1的shares
        for(int partyidx=0;partyidx<numofparties;++partyidx){
            //读取各个party的skN与invskN
            std::string party_dir = groupDir + "/party_" + std::to_string(partyidx);
            std::string skN_path = party_dir + "/skN";
            std::string invskN_path = party_dir + "/invskN";
            // 检查skN和invskN文件是否都存在
            if (!fs::exists(skN_path) || !fs::is_regular_file(skN_path) || 
                !fs::exists(invskN_path) || !fs::is_regular_file(invskN_path)) {
                std::cerr<<"party_"<<partyidx<<"的skN或invskN文件不存在"<<std::endl;
                break; 
            }
         
            std::ifstream skN_ifs(skN_path);
            if (!skN_ifs.is_open()) {
                std::cerr << "Error: 无法打开party_"<<partyidx<<"的skN文件 " << skN_path << std::endl;
                break;
            }
            NativeVector skNVec(N, Q);
            for (uint32_t i = 0; i < N; ++i) {
                std::string val_str;
                if (!(skN_ifs >> val_str)) {
                    break; 
                }
                skNVec[i] = NativeInteger(val_str);
            }
            skN_ifs.close();

            std::ifstream invskN_ifs(invskN_path);
            if (!invskN_ifs.is_open()) {
                std::cerr << "Error: 无法打开party_"<<partyidx<<"的invskN文件 " << invskN_path << std::endl;
                break;
            }
            NativeVector invskNVec(N, Q);
            for (uint32_t i = 0; i < N; ++i) {
                std::string val_str;
                if (!(invskN_ifs >> val_str)) {
                    break; 
                }
                invskNVec[i] = NativeInteger(val_str);
            }
            invskN_ifs.close();

            NativePoly skNPoly(cc.GetParams()->GetVectorNTRUParams()->GetPolyParams(), Format::COEFFICIENT, true);
            NativePoly invskNPoly(cc.GetParams()->GetVectorNTRUParams()->GetPolyParams(), Format::COEFFICIENT, true);
            skNPoly.SetValues(skNVec, Format::COEFFICIENT);
            invskNPoly.SetValues(invskNVec, Format::COEFFICIENT);
            skNPoly.SetFormat(Format::EVALUATION);      //注意，先转换成NTT形式再进行秘密共享
            invskNPoly.SetFormat(Format::EVALUATION);
        
            auto skN_share = secretSharing.Split(skNPoly, numofparties);
            auto invskN_share = secretSharing.Split(invskNPoly, numofparties);
            if(partyidx == 0){
                this->F_shares = skN_share;
                this->F_inv_shares = invskN_share;
            }else{
                this->F_shares = secretSharing.Mult(this->F_shares, skN_share, a_shares, b_shares, c_shares);
                this->F_inv_shares = secretSharing.Mult(this->F_inv_shares, invskN_share, a_shares, b_shares, c_shares);
            }
        }   

        //各参与方P_j的sk(s_j,0,s_j,1,...,s_j,n-1)，目标计算X^(-s_l)，s_l=\sum s_(j,l) 
        std::vector<NativeVector> sk_all_vec(numofparties, NativeVector(n, cc.GetParams()->GetLWEParams()->Getq()));
        for(int partyidx=0;partyidx<numofparties;++partyidx){
            //读入全部party的sk
            std::string party_sk_path = groupDir + "/party_" + std::to_string(partyidx) + "/sk";
            if (!fs::exists(party_sk_path)) {
                continue;
            }
            std::ifstream sk_ifs(party_sk_path);
            if (!sk_ifs.is_open()) {
                continue;
            }
            std::vector<std::string> values;
            std::string val_str;
            while (sk_ifs >> val_str) {
                values.push_back(val_str);
            }
            sk_ifs.close();
            if (values.empty()) {
                std::cerr<<"party_"<<partyidx<<"的sk文件为空"<<std::endl;
                break;
            }
            NativeVector sk_vec(n, Q);
            for (size_t i = 0; i < values.size(); ++i) {
                sk_all_vec[partyidx][i] = NativeInteger(values[i]);
            }
        }
        

        //计算第三个分片
        for(int i = 0; i<n; ++i){
            for(int partyidx=0;partyidx<numofparties;++partyidx){
                NativePoly zeroPoly(cc.GetParams()->GetVectorNTRUParams()->GetPolyParams(), Format::COEFFICIENT, true);
                if(sk_all_vec[partyidx][i] == 1){
                    zeroPoly[1] = 1;
                }else if(sk_all_vec[partyidx][i] == cc.GetParams()->GetLWEParams()->Getq()-1){
                    zeroPoly[N-1] = cc.GetParams()->GetLWEParams()->GetQ() - 1;
                }else if(sk_all_vec[partyidx][i] == 0){
                    zeroPoly[0] = 1;
                }else{  
                    cerr<<"party_"<<partyidx<<"的sk非法"<<endl;
                    break;
                }
                NativePoly Xpows(zeroPoly);
                Xpows.SetFormat(Format::EVALUATION);
                auto xpows_share = secretSharing.Split(Xpows, numofparties);
                if(partyidx == 0){
                    this->sk_sum_shares[i] = xpows_share;
                }else{
                    this->sk_sum_shares[i] = secretSharing.Mult(this->sk_sum_shares[i], xpows_share, a_shares, b_shares, c_shares);
                }
            }
        }

        //计算第四个分片
        for(int i = 0; i<n; ++i){
            for(int partyidx=0;partyidx<numofparties;++partyidx){
                NativePoly zeroPoly(cc.GetParams()->GetVectorNTRUParams()->GetPolyParams(), Format::COEFFICIENT, true);
                //注意取反
                if(sk_all_vec[partyidx][i] == 1){
                    zeroPoly[N-1] = cc.GetParams()->GetLWEParams()->GetQ() - 1;
                }else if(sk_all_vec[partyidx][i] == cc.GetParams()->GetLWEParams()->Getq()-1){
                    zeroPoly[1] = 1;
                }else if(sk_all_vec[partyidx][i] == 0){
                    zeroPoly[0] = 1;
                }else{  
                    cerr<<"party_"<<partyidx<<"的sk非法"<<endl;
                    break;
                }
                NativePoly Xpows(zeroPoly);
                Xpows.SetFormat(Format::EVALUATION);
                auto xpows_share = secretSharing.Split(Xpows, numofparties);
                if(i == 0 && partyidx == 0){
                    this->sk_inv_all_sum_shares = xpows_share;
                }else{
                    this->sk_inv_all_sum_shares = secretSharing.Mult(this->sk_inv_all_sum_shares, xpows_share, a_shares, b_shares, c_shares);
                }
            }
        }

    }
    
    //Getters
    const std::vector<lbcrypto::NativePoly>& GetF_shares() const { return F_shares; }
    const std::vector<lbcrypto::NativePoly>& GetF_inv_shares() const { return F_inv_shares; }
    const std::vector<std::vector<lbcrypto::NativePoly>>& Getsk_inv_sum_shares() const { return sk_sum_shares; }
    const std::vector<lbcrypto::NativePoly>& Getsk_all_sum_shares() const { return sk_inv_all_sum_shares; }
    
    

};

}

#endif // ENROLL_SHARES_H
