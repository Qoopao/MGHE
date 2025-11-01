#include "binfhecontext.h"
#include <string>
#include <vector>

#include <iostream>
#include <ctime>
#include <random>
 
using namespace lbcrypto;
using namespace std;

int main() {

    auto cc = BinFHEContext();
    cc.GenerateBinFHEContext(BIGP, XZDDF);

    //组数以及组内人数
    int k = 1;
    int ki = 1;
    //注意要先生成这个目录，否则会出问题
    cc.inituser(k, ki);

    std::vector<std::vector<LWEPrivateKey>> sk(k, std::vector<LWEPrivateKey>(ki));
    for(int i = 0; i < k; ++i){
        for(int j = 0; j < ki; ++j){
            auto lwesk = cc.LWESecretKeyGen(i, j);
            auto lwepk = cc.LWEPublicKeyGen(lwesk, i, j);
            sk[i][j]=lwesk;
        }
        auto jointpk = cc.LWEJointPublicKeyGen(cc.GetParams()->GetLWEParams(), i);
        
    }

    //多组加密
    //randomly generate a number from 0 to 1 
    std::random_device rd; 
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);
    for(int i=0;i<100;++i){
        LWEPlaintext m = dis(gen);
        LWECiphertext ct = cc.MGHE_Encrypt(m, 0);

        // 调试输出：加密后a、b
        // const auto& avec = ct->GetA();
        // std::cout << "加密: m=" << m << std::endl;
        // std::cout << "a = [";
        // for (size_t i = 0; i < avec.GetLength(); ++i) {
        //     std::cout << avec[i];
        //     if (i != avec.GetLength() - 1) std::cout << ", ";
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "b = " << ct->GetB() << std::endl;

        //多组解密
        LWEPlaintext result;
        //cc.MGHE_Decrypt(sk, ct, &result);=下面这个大括号的内容
        {
            auto params = cc.GetParams()->GetLWEParams();
            LWEPlaintextModulus p = 4;
            const NativeInteger& mod = ct->GetModulus();
            if (mod % (p * 2) != 0 && mod.ConvertToInt() & (1 == 0)) {
                std::string errMsg = "ERROR: ciphertext modulus q needs to be divisible by plaintext modulus p*2.";
                OPENFHE_THROW(not_implemented_error, errMsg);
            }
        
            NativeVector a   = ct->GetA();
            int k=sk.size();
            int ki=sk[0].size();
            // std::cout<<"DECRYPT: a = "<<a<<std::endl;
            // std::cout<<"DECRYPT: k = "<<k<<", ki = "<<ki<<std::endl;
        
            
            std::vector<std::vector<NativeVector>> s(k,std::vector<NativeVector>(ki));
            for(int i=0;i<k;++i){
                for(int j=0;j<ki;++j){
                    s[i][j]=sk[i][j]->GetElement();
                }
            }
            //std::cout<<"DECRYPT: s[0][0] = "<<s[0][0]<<std::endl;
        
        
            NativeInteger mu = mod.ComputeMu();
            for(int i=0;i<k;++i){
                for(int j=0;j<ki;++j){
                    s[i][j].SwitchModulus(mod);
                }
            }
        
            NativeInteger inner(0);
            for(int i = 0; i < k; ++i){
                for (int j = 0; j < ki; ++j) {
                    for(uint32_t l=0;l<params->Getn();++l){
                        inner += a[params->Getn()*i+l].ModMulFast(s[i][j][l], mod, mu);
                    }
                }
            }
            //std::cout<<"DECRYPT: inner before ModEq = "<<inner<<std::endl;
        
            inner.ModEq(mod);                                       //=inner%mod
            //std::cout<<"DECRYPT: inner after ModEq = "<<inner<<std::endl;
            NativeInteger r = ct->GetB();   
            //std::cout<<"DECRYPT: b (ct->GetB()) = "<<r<<std::endl;
            r.ModAddFastEq(inner, mod);                             // r = B + \sum{ai*si} % q
            //std::cout<<"DECRYPT: r after b + inner = "<<r<<std::endl;
            // 使用标准LWE解密的舍入方法，正确处理负数情况
            //std::cout<<"DECRYPT: final r before rounding = "<<r<<std::endl;
            //std::cout<<"DECRYPT: 2*p = "<<(2*p)<<std::endl;
            
            // // 添加 q/(2p) 进行舍入调整
            // r.ModAddFastEq((mod / (p * 2)), mod);
            // std::cout<<"DECRYPT: r after adding q/(2p) = "<<r<<std::endl;
            
            // // 使用整数除法进行舍入
            // *result = ((NativeInteger(p) * r) / mod).ConvertToInt();
            // std::cout<<"DECRYPT: final result = "<<*result<<std::endl;
        
            if(r>mod/2){
                r=mod-r;
            }
            //std::cout<<"DECRYPT: final result = "<<*result<<std::endl;
        
            result = (r.MultiplyAndRound(NativeInteger(4),mod)).ConvertToInt();
        }



        std::cout << "result=" << result << ", m=" << m;
        if(result==m){
            std::cout << " ✅" << std::endl;
        }else{
            std::cout << " ❌" << std::endl;
        }
    }


    return 0;
}