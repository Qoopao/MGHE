#include <bits/stdint-uintn.h>
#include "binfhe-constants.h"
#include "binfhecontext.h"
#include "lattice/hal/lat-backend.h"
#include "lwe-ciphertext-fwd.h"
#include "lwe-privatekey.h"
#include "lwe-publickey-fwd.h"
#include "math/hal/nativeintbackend.h"
#include "utils/inttypes.h"
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <system_error>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


using namespace lbcrypto;
namespace fs = std::filesystem;

//只需加一个参数（密文模数），因为LWE维度不变
std::vector<NativeVector> LWECRMGen(BinFHEContext& cc, NativeInteger modulus){
    
    auto LWEParams = cc.GetParams()->GetLWEParams();
    auto n = LWEParams->Getn();
    
    // 生成CRM矩阵
    std::vector<NativeVector> crm;
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(modulus);
    for (usint i = 0; i < n; ++i) {
        crm.push_back(dug.GenerateVector(n, modulus));
    }

    return crm;

}

LWEPrivateKey GENLWESK(BinFHEContext& cc){
    auto LWEParams = cc.GetParams()->GetLWEParams();
    auto modulus = cc.GetParams()->GetLWEParams()->Getq();
    auto n = LWEParams->Getn();
    TernaryUniformGeneratorImpl<NativeVector> tug;
    
    NativeVector sk_vec = tug.GenerateVector(n, modulus);
    LWEPrivateKey sk = std::make_shared<LWEPrivateKeyImpl>(sk_vec);
    return sk;
}

LWEPublicKey GENLWEPK(const std::shared_ptr<LWECryptoParams>& params, std::vector<NativeVector> CRM, ConstLWEPrivateKey& sk, NativeInteger modulus)  {
    size_t dim = params->Getn();
    
    // Generate error vector e
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    NativeVector e = dgg.GenerateVector(dim, modulus);
    //std::cout<<"e="<<e<<std::endl;

    // Compute r = -Ms + e
    NativeVector r(dim, modulus);
    NativeVector ske = sk->GetElement();
    ske.SwitchModulus(modulus);
    //cout<<"ske="<<ske<<endl;
    NativeInteger mu = modulus.ComputeMu();
    for (size_t j = 0; j < dim; ++j) {
        NativeInteger sum = 0;
        for (size_t i = 0; i < dim; ++i) {
            sum.ModAddFastEq(CRM[j][i].ModMulFast(ske[i], modulus, mu), modulus);
        }
        // r[j] = (e[j] + (modulus - sum)) % modulus
        r[j] = (e[j] + (modulus - sum.Mod(modulus))).Mod(modulus);
    }

    return std::make_shared<LWEPublicKeyImpl>(LWEPublicKeyImpl(std::move(CRM), std::move(r)));
}

LWECiphertext ENCLWE(BinFHEContext& cc, std::vector<NativeVector> M, LWEPublicKey PK, int m){
    
    auto LWEParams = cc.GetParams()->GetLWEParams();
    size_t n = LWEParams->Getn();
    auto modulus = M[0].GetModulus();

    // 3. 生成随机向量v ∈ Zq^n (从二进制分布)
    BinaryUniformGeneratorImpl<NativeVector> bug;
    NativeVector v = bug.GenerateVector(n, modulus);
    //std::cout<<"v="<<v<<std::endl;

    // 4. 采样噪声e1, e2
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    double STD_DEV = 3;
    dgg.SetStd(STD_DEV);
    NativeInteger e1 = dgg.GenerateInteger(modulus); // 标量噪声
    NativeVector e2 = dgg.GenerateVector(n, modulus); // 向量噪声
    //cout<<"e1="<<e1<<endl;
    //std::cout<<"e2="<<e2<<std::endl;
    

    // 5. 计算b = <jointpk, v> + (q/4)*m + e1
    NativeInteger b = 0;
    
    NativeInteger mu = modulus.ComputeMu();
    for (size_t i = 0; i < n; ++i) {
        b.ModAddFastEq(PK->Getv()[i].ModMulFast(v[i], modulus, mu), modulus);
    }
    //std::cout<<"b after <pk,v> = "<<b<<std::endl;

    b.ModAddFastEq((m % 4) * (modulus / 4), modulus);
    
    
    b.ModAddFastEq(e1, modulus);
    //std::cout<<"b after adding e1 = "<<b<<std::endl;
    // 6. 计算a = v^T * M + e2
    NativeVector a(n, modulus);
    for (size_t j = 0; j < n; ++j) {
        NativeInteger sum = 0;
        for (size_t i = 0; i < n; ++i) {
            sum.ModAddFastEq(v[i].ModMulFast(M[i][j], modulus, mu), modulus);
        }
        a[j] = sum.ModAdd(e2[j], modulus);
    }
    
    // 7. 输出密文(a, b)
    auto ct = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b.Mod(modulus)));
    
    ct->SetptModulus(4);
    return ct;
}

// 修改DECLWE函数，返回e的绝对值供统计使用
LWEPlaintext DECLWE(BinFHEContext& cc, LWECiphertext ct, LWEPrivateKey sk, int64_t& e_abs){
    auto LWEParams = cc.GetParams()->GetLWEParams();
    size_t n = LWEParams->Getn();
    auto modulus = LWEParams->Getq();

    NativeVector a   = ct->GetA();
    NativeVector s = sk->GetElement();
    s.SwitchModulus(modulus);
    //cout<<"s="<<s<<endl;

    NativeInteger mu = modulus.ComputeMu();
    NativeInteger inner(0);
    for(uint32_t l=0;l<n;++l){
        inner += a[l].ModMulFast(s[l], modulus, mu);
    }
    inner.ModEq(modulus);                                        
    
    NativeInteger r = ct->GetB();   
    
    r.ModAddFastEq(inner, modulus);                            

    auto rcopy = r;
    //cout<<"rcopy="<<rcopy<<" "<<endl;

    if(r>modulus/2){
        r = modulus-r;
    }
    auto res = (r.MultiplyAndRound(NativeInteger(4),modulus));
    LWEPlaintext result = res.ConvertToInt();

    NativeInteger e;
    if(res==1){
        //b+as = (q/4) + e
        // 处理无符号整数减法可能的下溢问题
        if(rcopy >= modulus/4) {
            e = rcopy - (modulus/4);
        } else {
            // 当rcopy < modulus/4时，e应该是负数
            // 但由于NativeInteger是无符号的，我们需要特殊处理
            // 这里先记录为小数值，后续处理绝对值时再转换
            e = rcopy;
        }
        //cout<<"res==1时的e="<<e<<endl;
    }else{
        //b+as = e
        e = rcopy;
    }
    // if(rcopy > modulus/2){
    //     cout<<"解密时r大于q/2, ropy="<<rcopy<<", e="<<e<<endl;
    // }
    
    // 计算e的绝对值
    if(res == 1) {
        // 对于res==1的情况，需要特殊处理e的绝对值
        if(rcopy >= modulus/4) {
            // 正数情况：e = rcopy - modulus/4
            e_abs = e.ConvertToInt();
            // 打印大于39的e值
            //if(e_abs > 39){
                cout<<"e="<<e_abs<<endl;
            //}
        } else {
            // 负数情况：e = rcopy - modulus/4 应该是负数
            // 计算其绝对值：modulus/4 - rcopy
            e_abs = ((modulus/4) - rcopy).ConvertToInt();
            // 打印大于39的e值
            //if(e_abs > 39){
                cout<<"e=-"<<e_abs<<endl;
            //}
        }
    } else {
        // 对于res==0的情况，使用原来的逻辑
        if( e > modulus/2){
            //cout<<"modulus/4= "<<modulus/4<<" ";
            //cout<<"ori e="<<e<<endl;
            e_abs = (modulus - e).ConvertToInt();
            // 打印大于39的e值
            if(e_abs > 39){
                cout<<"e=-"<<e_abs<<endl;
            }
        }else{
            e_abs = e.ConvertToInt();
            // 打印大于39的e值
            //if(e_abs > 39){
                cout<<"e="<<e_abs<<endl;
            //}
        }
    }

    return result;
}
 
    

int main(){
    NativeInteger q1 = FreshEncModulus;
    auto cc = BinFHEContext();
    cc.GenerateBinFHEContext(BIGP, XZDDF);

    auto crm = LWECRMGen(cc, q1);
    auto sk = GENLWESK(cc);
    auto pk = GENLWEPK(cc.GetParams()->GetLWEParams(), crm, sk, q1);

    // 计数器初始化
    int total_tests = 10;  // 运行一万次测试
    int fail_count = 0;
    int large_error_count = 0;
    const int64_t ERROR_THRESHOLD = 39;

    cout << "开始运行" << total_tests << "次测试..." << endl;
    cout << "--------------------------------------------------" << endl;

    for(int i=0;i<total_tests;++i){
        // 随机生成0或1
        int m = rand() % 2;  // 也可以改为 rand() % 2 来随机生成0或1
        auto ct = ENCLWE(cc, crm, pk, m);
        auto ctMS = cc.GetLWEScheme()->ModSwitch(cc.GetParams()->GetLWEParams()->Getq(), ct);
        
        int64_t e_abs = 0;
        auto res = DECLWE(cc, ctMS, sk, e_abs);
        
        // 统计错误
        if(res != m){
            //cout << "fail" << endl;
            fail_count++;
        }else{
            //cout << "success" << endl;
        }
        

        // 统计e绝对值大于39的情况
        if(e_abs > ERROR_THRESHOLD){
            large_error_count++;
        }
        
        // 每1000次输出一次进度
        if((i + 1) % 10 == 0){
            cout << "已完成 " << (i + 1) << " 次测试..." << endl;
        }
    }
    
    cout << "--------------------------------------------------" << endl;
    cout << "测试完成！" << endl;
    cout << "总测试次数: " << total_tests << endl;
    cout << "失败次数: " << fail_count << endl;
    cout << "失败概率: " << (double)(fail_count * 100 / total_tests) << "%" << endl;
    cout << "e绝对值>" << ERROR_THRESHOLD << "的次数: " << large_error_count << endl;
    cout << "e绝对值>" << ERROR_THRESHOLD << "的概率: " << (double)(large_error_count * 100 / total_tests) << "%" << endl;
    
    return 0;
}