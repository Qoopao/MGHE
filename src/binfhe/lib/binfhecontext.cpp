#include "binfhecontext.h"
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <filesystem>
#include <fstream>
#include "cereal/archives/json.hpp"
#include "lattice/hal/lat-backend.h"
#include "lwe-ciphertext-fwd.h"
#include "math/hal/nativeintbackend.h"
#include "utils/inttypes.h"
#include "rlwe-privatekey-fwd-lyh.h"
#include "rlwe-privatekey-lyh.h"
#include "rlwe-publickey-fwd-lyh.h"
#include "rlwe-publickey-lyh.h"

namespace lbcrypto {

void BinFHEContext::GenerateBinFHEContext(uint32_t n, uint32_t N, const NativeInteger& q, const NativeInteger& Q,
                                          double std, uint32_t baseKS, uint32_t baseG, uint32_t baseR,
                                          SecretKeyDist keyDist, BINFHE_METHOD method, uint32_t numAutoKeys) {
    //用 std::make_shared 创建了 std::shared_ptr 的智能指针实例
    auto lweparams = std::make_shared<LWECryptoParams>(n, N, q, Q, Q, std, baseKS);
    auto rgswparams =
        std::make_shared<RingGSWCryptoParams>(N, Q, q, baseG, baseR, method, std, keyDist, true, numAutoKeys);
    m_params       = std::make_shared<BinFHECryptoParams>(lweparams, rgswparams);
    m_binfhescheme = std::make_shared<BinFHEScheme>(method);
}

void BinFHEContext::GenerateBinFHEContext(BINFHE_PARAMSET set, bool arbFunc, uint32_t logQ, int64_t N,
                                          BINFHE_METHOD method, bool timeOptimization) {
    if (GINX != method) {
        std::string errMsg("ERROR: CGGI is the only supported method");
        OPENFHE_THROW(not_implemented_error, errMsg);
    }
    if (set != STD128 && set != TOY) {
        std::string errMsg("ERROR: STD128 and TOY are the only supported sets");
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    if (logQ > 29) {
        std::string errMsg("ERROR: logQ > 29 is not supported");
        OPENFHE_THROW(not_implemented_error, errMsg);
    }
    if (logQ < 11) {
        std::string errMsg("ERROR: logQ < 11 is not supported");
        OPENFHE_THROW(not_implemented_error, errMsg);
    }
    auto logQprime = 54;
    uint32_t baseG = 0;
    if (logQ > 25) {
        baseG = 1 << 14;
    }
    else if (logQ > 16) {
        baseG = 1 << 18;
    }
    else if (logQ > 11) {
        baseG = 1 << 27;
    }
    else {  // if (logQ == 11)
        baseG     = 1 << 5;
        logQprime = 27;
    }

    m_timeOptimization = timeOptimization;
    SecurityLevel sl   = HEStd_128_classic;
    // choose minimum ringD satisfying sl and Q
    uint32_t ringDim = StdLatticeParm::FindRingDim(HEStd_ternary, sl, logQprime);
    if (N >= ringDim) {  // if specified some larger N, security is also satisfied
        ringDim = N;
    }
    // find prime Q for NTT
    NativeInteger Q = PreviousPrime<NativeInteger>(FirstPrime<NativeInteger>(logQprime, 2 * ringDim), 2 * ringDim);
    // q = 2*ringDim by default for maximum plaintext space, if needed for arbitrary function evaluation, q = ringDim
    uint32_t q = arbFunc ? ringDim : 2 * ringDim;

    uint64_t qKS = 1 << 30;
    qKS <<= 5;

    uint32_t n      = (set == TOY) ? 32 : 1305;
    auto lweparams  = std::make_shared<LWECryptoParams>(n, ringDim, q, Q, qKS, 3.19, 32);
    auto rgswparams = std::make_shared<RingGSWCryptoParams>(ringDim, Q, q, baseG, 23, method, 3.19, UNIFORM_TERNARY,
                                                            ((logQ != 11) && timeOptimization));

    m_params       = std::make_shared<BinFHECryptoParams>(lweparams, rgswparams);
    m_binfhescheme = std::make_shared<BinFHEScheme>(method);

#if defined(BINFHE_DEBUG)
    std::cout << ringDim << " " << Q < < < < " " << n << " " << q << " " << baseG << std::endl;
#endif
}

//用这个
void BinFHEContext::GenerateBinFHEContext(BINFHE_PARAMSET set, BINFHE_METHOD method) {
    enum { PRIME = 0 };  // value for modKS if you want to use the intermediate prime for modulus for key switching
    constexpr double STD_DEV = 3 ;
    constexpr double STD_NTRU = 0.5;
    // clang-format off

    /*
    gadgetBase就是baseG
    */
    const std::unordered_map<BINFHE_PARAMSET, BinFHEContextParams> paramsMap
    ({
        //               numberBits|cyclOrder|latticeParam|  mod|   modKS|  stdDev| baseKS| gadgetBase| baseRK| numAutoKeys| keyDist
        { TOY,               { 27,     1024,          64,  512,   PRIME, STD_DEV,     25,    1 <<  9,  23,     9,  UNIFORM_TERNARY} },
        { MEDIUM,            { 28,     2048,         422, 1024, 1 << 14, STD_DEV, 1 << 7,    1 << 10,  32,    10,  UNIFORM_TERNARY} },
        { STD256,            { 29,     4096,         990, 2048, 1 << 14, STD_DEV, 1 << 7,    1 <<  8,  46,    10,  UNIFORM_TERNARY} },
        { STD128Q,           { 25,     2048,         534, 1024, 1 << 14, STD_DEV,     32,    1 <<  7,  32,    10,  UNIFORM_TERNARY} },
        { STD128Q_LMKCDEY,   { 27,     2048,         448, 1024, 1 << 13, STD_DEV,     32,    1 <<  9,  32,    10,  GAUSSIAN       } },
        { STD192Q,           { 35,     4096,         875, 1024, 1 << 15, STD_DEV,     32,    1 << 12,  32,    10,  UNIFORM_TERNARY} },
        { STD256Q,           { 27,     4096,        1225, 1024, 1 << 16, STD_DEV,     16,    1 <<  7,  32,    10,  UNIFORM_TERNARY} },
        { STD128_3,          { 27,     2048,         541, 1024, 1 << 15, STD_DEV,     32,    1 <<  7,  32,    10,  UNIFORM_TERNARY} },
        { STD128_3_LMKCDEY,  { 28,     2048,         485, 1024, 1 << 15, STD_DEV,     32,    1 << 10,  32,    10,  GAUSSIAN       } },
        { STD128Q_3,         { 50,     4096,         575, 2048, 1 << 15, STD_DEV,     32,    1 << 25,  32,    10,  UNIFORM_TERNARY} },
        { STD128Q_3_LMKCDEY, { 27,     2048,         524, 1024, 1 << 15, STD_DEV,     32,    1 <<  9,  32,    10,  GAUSSIAN       } },
        { STD192Q_3,         { 34,     4096,         922, 2048, 1 << 16, STD_DEV,     16,    1 << 12,  32,    10,  UNIFORM_TERNARY} },
        { STD256Q_3,         { 27,     4096,        1400, 4096, 1 << 16, STD_DEV,     21,    1 <<  6,  32,    10,  UNIFORM_TERNARY} },
        { STD128_4,          { 27,     2048,         541, 2048, 1 << 15, STD_DEV,     32,    1 <<  7,  32,    10,  UNIFORM_TERNARY} },
        { STD128_4_LMKCDEY,  { 28,     2048,         522, 2048, 1 << 15, STD_DEV,     32,    1 << 10,  32,    10,  GAUSSIAN       } },
        { STD128Q_4,         { 50,     4096,         647, 2048, 1 << 16, STD_DEV,     16,    1 << 25,  32,    10,  UNIFORM_TERNARY} },
        { STD128Q_4_LMKCDEY, { 27,     2048,         524, 2048, 1 << 15, STD_DEV,     32,    1 <<  7,  32,    10,  GAUSSIAN       } },
        { STD192Q_4,         { 34,     4096,         980, 2048, 1 << 17, STD_DEV,     16,    1 << 12,  32,    10,  UNIFORM_TERNARY} },
        { STD256Q_4,         { 27,     4096,        1625, 4096, 1 << 21, STD_DEV,     16,    1 <<  6,  32,    10,  UNIFORM_TERNARY} },
        { SIGNED_MOD_TEST,   { 28,     2048,         512, 1024,   PRIME, STD_DEV,     25,    1 <<  7,  23,    10,  UNIFORM_TERNARY} },
        //               numberBits|cyclOrder|latticeParam|  mod|   modKS|  stdDev| baseKS| gadgetBase| baseRK| numAutoKeys| keyDist
        { STD128_LMKCDEY,    { 28,     2048,         446, 1024, 1 << 13, STD_DEV, 1 << 5,    1 << 10,  32,    10,  GAUSSIAN       } },
        { STD128_LMKCDEY_New,{ 28,     2048,         446, 1024, 1 << 13, STD_DEV, 1 << 5,    1 <<  7,  32,    10,  GAUSSIAN       } },
        { STD128_AP,         { 27,     2048,         503, 1024, 1 << 14, STD_DEV, 1 << 5,    1 <<  9,  32,    10,  UNIFORM_TERNARY} },
        { STD128,            { 27,     2048,         503, 1024, 1 << 14, STD_DEV, 1 << 5,    1 <<  9,  32,    10,  UNIFORM_TERNARY} },
        //               numberBits|cyclOrder|latticeParam|  mod|   modKS|  stdDev| baseKS| gadgetBase| baseRK| numAutoKeys| keyDist
        { P128T,             { 21,     2048,         512, 1024, 1 << 14, STD_NTRU,    32,    1 <<  7,  32,    10,  UNIFORM_TERNARY} },
        { P128G,             { 21,     2048,         446, 1024, 1 << 14, STD_NTRU,    32,    1 <<  7,  32,    10,  GAUSSIAN       } },
        { P128T_2,           { 21,     2048,         512, 1024, 1 << 14, STD_NTRU,    32,    1 <<  6,  32,    10,  UNIFORM_TERNARY} },
        { P128G_2,           { 21,     2048,         446, 1024, 1 << 14, STD_NTRU,    32,    1 <<  6,  32,    10,  GAUSSIAN       } },
        //               numberBits|cyclOrder|latticeParam|  mod|   modKS|  stdDev| baseKS| gadgetBase| baseRK| numAutoKeys| keyDist
        { STD192,            { 37,     4096,         805, 1024, 1 << 15, STD_DEV,     32,    1 << 13,  32,    10,  UNIFORM_TERNARY} },
        //               numberBits|cyclOrder|latticeParam|  mod|   modKS|  stdDev| baseKS| gadgetBase| baseRK| numAutoKeys| keyDist
        { P192T,             { 26,     4096,        1024, 1024, 1 << 17, STD_NTRU,    28,    1 <<  9,  32,    10,  UNIFORM_TERNARY} },
        { P192G,             { 26,     4096,         805, 1024, 1 << 17, STD_NTRU,    28,    1 <<  9,  32,    10,  GAUSSIAN       } },
        //                               Q,                 cyclOrder,                  n,             q,    q_ks,   std_dev,  B_ks,    B
        { BIGP,              { 60,     4096,        16, 2048, 1 << 18, STD_NTRU,   32,    1 <<  10,  32,    10,  UNIFORM_TERNARY} },
        
        { L100T_2,           { 50,     8192,        600,  4096, 1 << 30, STD_NTRU,   1024,    1 <<  5,  32,    10,  UNIFORM_TERNARY} },
        { L100T_4,           { 50,     8192,        600,  4096, 1 << 30, STD_NTRU,   1024,    1 <<  5,  32,    10,  UNIFORM_TERNARY} },
        { L100T_8,           { 50,     8192,        600, 4096, 1 << 30, STD_NTRU,   1024,    1 <<  4,  32,    10,  UNIFORM_TERNARY} },

        { L128T_2,           { 60,     16384,        600,  4096, 1 << 30, STD_NTRU,   1024,    1 <<  9,  32,    10,  UNIFORM_TERNARY} },
        { L128T_4,           { 60,     16384,        600,  4096, 1 << 30, STD_NTRU,   1024,    1 <<  9,  32,    10,  UNIFORM_TERNARY} },
        { L128T_8,           { 60,     16384,        600, 4096, 1 << 30, STD_NTRU,   1024,    1 <<  8,  32,    10,  UNIFORM_TERNARY} },
    });//BIGP,              { 50,     1024,         512, 1024, 1 << 13, STD_NTRU,   28,    1 <<  9,  32,    10,  }
    // clang-format on

    auto search = paramsMap.find(set);
    if (paramsMap.end() == search) {
        std::string errMsg("ERROR: Unknown parameter set [" + std::to_string(set) + "] for FHEW.");
        OPENFHE_THROW(config_error, errMsg);
    }

    BinFHEContextParams params = search->second;
    // intermediate prime


    // 1. 首先计算第一个参数：获取具有指定比特数的第一个素数
    NativeInteger firstPrime = FirstPrime<NativeInteger>(params.numberBits, params.cyclOrder);
    // 2. 然后计算第二个参数：获取比firstPrime小的前一个素数
    NativeInteger previousPrime = PreviousPrime<NativeInteger>(firstPrime, params.cyclOrder);
    // 3. 最后用这两个值构造Q对象，最终NTRU/RLWE密文模数用的这个
    NativeInteger Q(previousPrime);

    usint ringDim  = params.cyclOrder / 2;

    //这里设置KS时是否使用一个中间模数
    auto lweparams = (PRIME == params.modKS) ?
                         std::make_shared<LWECryptoParams>(params.latticeParam, ringDim, params.mod, Q, Q,
                                                           params.stdDev, params.baseKS, params.keyDist) :
                         std::make_shared<LWECryptoParams>(params.latticeParam, ringDim, params.mod, Q, params.modKS,
                                                           params.stdDev, params.baseKS, params.keyDist);
    
    
    
    
    if(method == XZDDF)
    {
        auto vntruparams =
        std::make_shared<VectorNTRUCryptoParams>(ringDim, Q, params.mod, params.gadgetBase, params.baseRK, method,
                                              params.stdDev, params.keyDist, false, params.numAutoKeys);
        m_params       = std::make_shared<BinFHECryptoParams>(lweparams, vntruparams);
        
    }else{
        auto rgswparams =
        std::make_shared<RingGSWCryptoParams>(ringDim, Q, params.mod, params.gadgetBase, params.baseRK, method,
                                              params.stdDev, params.keyDist, false, params.numAutoKeys);
        m_params       = std::make_shared<BinFHECryptoParams>(lweparams, rgswparams);
    }//wkx

    
    m_binfhescheme = std::make_shared<BinFHEScheme>(method);

}

void BinFHEContext::GenerateBinFHEContext(const BinFHEContextParams& params, BINFHE_METHOD method) {
    enum { PRIME = 0 };  // value for modKS if you want to use the intermediate prime for modulus for key switching
    // intermediate prime
    NativeInteger Q(
        PreviousPrime<NativeInteger>(FirstPrime<NativeInteger>(params.numberBits, params.cyclOrder), params.cyclOrder));

    usint ringDim = params.cyclOrder / 2;

    auto lweparams = (PRIME == params.modKS) ?
                         std::make_shared<LWECryptoParams>(params.latticeParam, ringDim, params.mod, Q, Q,
                                                           params.stdDev, params.baseKS, params.keyDist) :
                         std::make_shared<LWECryptoParams>(params.latticeParam, ringDim, params.mod, Q, params.modKS,
                                                           params.stdDev, params.baseKS, params.keyDist);

    auto rgswparams =
        std::make_shared<RingGSWCryptoParams>(ringDim, Q, params.mod, params.gadgetBase, params.baseRK, method,
                                              params.stdDev, params.keyDist, false, params.numAutoKeys);

    m_params       = std::make_shared<BinFHECryptoParams>(lweparams, rgswparams);
    m_binfhescheme = std::make_shared<BinFHEScheme>(method);
}

LWEPrivateKey BinFHEContext::KeyGen() const {
    auto& LWEParams = m_params->GetLWEParams(); //m_params是指向BinFHECryptoParams对象的共享指针,GetLWEParams()返回的是指向LWECryptoParams对象的共享指针
    //根据LWE参数中"分布"这个枚举变量，选择合适的分布
    if (LWEParams->GetKeyDist() == GAUSSIAN)
        return m_LWEscheme->KeyGenGaussian(LWEParams->Getn(), LWEParams->GetqKS());
    return m_LWEscheme->KeyGen(LWEParams->Getn(), LWEParams->GetqKS());
}

LWEPrivateKey BinFHEContext::KeyGenN() const {
    auto& LWEParams = m_params->GetLWEParams();
    if (LWEParams->GetKeyDist() == GAUSSIAN)
        return m_LWEscheme->KeyGenGaussian(LWEParams->GetN(), LWEParams->GetQ());
    return m_LWEscheme->KeyGen(LWEParams->GetN(), LWEParams->GetQ());
}

LWEKeyPair BinFHEContext::KeyGenPair() const {
    auto&& LWEParams = m_params->GetLWEParams();
    return m_LWEscheme->KeyGenPair(LWEParams);
}

LWEPublicKey BinFHEContext::PubKeyGen(ConstLWEPrivateKey& sk) const {
    auto&& LWEParams = m_params->GetLWEParams();
    return m_LWEscheme->PubKeyGen(LWEParams, sk);
}

LWECiphertext BinFHEContext::Encrypt(ConstLWEPrivateKey& sk, LWEPlaintext m, BINFHE_OUTPUT output,
                                     LWEPlaintextModulus p, const NativeInteger& mod) const {
    const auto& LWEParams = m_params->GetLWEParams();

    LWECiphertext ct = (mod == 0) ? m_LWEscheme->Encrypt(LWEParams, sk, m, p, LWEParams->Getq()) :
                                    m_LWEscheme->Encrypt(LWEParams, sk, m, p, mod);

    // BINFHE_OUTPUT is kept as it is for backward compatibility but
    // this logic is obsolete now and commented out
    // if ((output != FRESH) && (p == 4)) {
    //    ct = m_binfhescheme->Bootstrap(m_params, m_BTKey, ct);
    //}

    return ct;
}



//多组加密
LWECiphertext BinFHEContext::MGHE_Encrypt(LWEPlaintext m, int groupidx, BINFHE_OUTPUT output,
                                LWEPlaintextModulus p, const NativeInteger& mod) const {
    auto&& LWEParams = m_params->GetLWEParams();
    return m_LWEscheme->MGHE_Encrypt(LWEParams, m, groupidx, p, LWEParams->Getq());
}

//测试用多组加密
LWECiphertext BinFHEContext::MGHE_ALLEncrypt(LWEPlaintext m, int k, BINFHE_OUTPUT output,
                                LWEPlaintextModulus p, const NativeInteger& mod) const {
    auto&& LWEParams = m_params->GetLWEParams();
    return m_LWEscheme->MGHE_ALLEncrypt(LWEParams, k, m, p, mod);
}

//多组解密
void BinFHEContext::MGHE_Decrypt(std::vector<std::vector<LWEPrivateKey>>& sk, ConstLWECiphertext& ct, LWEPlaintext* result, LWEPlaintextModulus p) const {
    auto&& LWEParams = m_params->GetLWEParams();
    m_LWEscheme->MGHE_Decrypt(LWEParams, sk, ct, result, p);
}










LWECiphertext BinFHEContext::Encrypt(ConstLWEPublicKey& pk, LWEPlaintext m, BINFHE_OUTPUT output, LWEPlaintextModulus p,
                                     const NativeInteger& mod) const {
    const auto& LWEParams = m_params->GetLWEParams();

    LWECiphertext ct = (mod == 0) ? m_LWEscheme->EncryptN(LWEParams, pk, m, p, LWEParams->GetQ()) :
                                    m_LWEscheme->EncryptN(LWEParams, pk, m, p, mod);

    // Switch from ct of modulus Q and dimension N to smaller q and n
    // This is done by default while calling Encrypt but the output could
    // be set to LARGE_DIM to skip this switching
    if (output == SMALL_DIM) {
        LWECiphertext ct1 = SwitchCTtoqn(m_BTKey.KSkey, ct);
        return ct1;
    }
    return ct;
}

LWECiphertext BinFHEContext::SwitchCTtoqn(ConstLWESwitchingKey& ksk, ConstLWECiphertext& ct) const {
    const auto& LWEParams = m_params->GetLWEParams();
    auto Q                = LWEParams->GetQ();
    auto N                = LWEParams->GetN();

    if ((ct->GetLength() != N) && (ct->GetModulus() != Q)) {
        std::string errMsg("ERROR: Ciphertext dimension and modulus are not large N and Q");
        OPENFHE_THROW(config_error, errMsg);
    }

    LWECiphertext ct1 = m_LWEscheme->SwitchCTtoqn(LWEParams, ksk, ct);

    return ct1;
}

void BinFHEContext::Decrypt(ConstLWEPrivateKey& sk, ConstLWECiphertext& ct, LWEPlaintext* result,
                            LWEPlaintextModulus p) const {
    auto&& LWEParams = m_params->GetLWEParams();
    m_LWEscheme->Decrypt(LWEParams, sk, ct, result, p);
}

LWESwitchingKey BinFHEContext::KeySwitchGen(ConstLWEPrivateKey& sk, ConstLWEPrivateKey& skN) const {
    return m_LWEscheme->KeySwitchGen(m_params->GetLWEParams(), sk, skN);
}

void BinFHEContext::BTKeyGen(ConstLWEPrivateKey& sk, KEYGEN_MODE keygenMode) {
    auto& RGSWParams = m_params->GetRingGSWParams();

    auto temp = RGSWParams->GetBaseG();

    //预先生成一些密钥
    if (m_timeOptimization) {
        auto gpowermap = RGSWParams->GetGPowerMap();
        for (std::map<uint32_t, std::vector<NativeInteger>>::iterator it = gpowermap.begin(); it != gpowermap.end();
             ++it) {
            RGSWParams->Change_BaseG(it->first);
            m_BTKey_map[it->first] = m_binfhescheme->KeyGen(m_params, sk, keygenMode);
        }
        RGSWParams->Change_BaseG(temp);
    }

    if (m_BTKey_map.size() != 0) {
        m_BTKey = m_BTKey_map[temp];
    }
    else {
        m_BTKey           = m_binfhescheme->KeyGen(m_params, sk, keygenMode);
        m_BTKey_map[temp] = m_BTKey;
    }
}

LWECiphertext BinFHEContext::EvalBinGate(const BINGATE gate, ConstLWECiphertext& ct1, ConstLWECiphertext& ct2) const {
    //std::cout<<"LWECiphertext BinFHEContext::EvalBinGate"<<std::endl;
    auto& VNTRUParams = m_params->GetVectorNTRUParams();

    if(VNTRUParams != nullptr)
    {   //需要改成m_NBTKey
        //std::cout<<"EvalBinGate XZDDF"<<std::endl;
        return m_binfhescheme->EvalBinGate(m_params, gate, m_NBTKey, ct1, ct2);
    }else{
        //std::cout<<"EvalBinGate AP or GINX"<<std::endl;
        return m_binfhescheme->EvalBinGate(m_params, gate, m_BTKey, ct1, ct2);
    }
}


LWECiphertext BinFHEContext::MGEvalBinGate(const BINGATE gate, ConstLWECiphertext& ct1, ConstLWECiphertext& ct2, int numofgroups) const {
    //std::cout<<"LWECiphertext BinFHEContext::EvalBinGate"<<std::endl;
    auto& VNTRUParams = m_params->GetVectorNTRUParams();
    if(VNTRUParams != nullptr)
    {   //需要改成m_NBTKey
        //std::cout<<"EvalBinGate XZDDF"<<std::endl;
        return m_binfhescheme->MGEvalBinGate(m_params, gate, m_NBTKey, ct1, ct2, numofgroups);
    }else{
        //std::cout<<"EvalBinGate AP or GINX"<<std::endl;
        //return m_binfhescheme->MGEvalBinGate(m_params, gate, m_BTKey, ct1, ct2);
        std::cout<<"VNTRUParams is empty"<<std::endl;
        return nullptr;
    }
}



LWECiphertext BinFHEContext::EvalBinGate(const BINGATE gate, const std::vector<LWECiphertext>& ctvector) const {
    return m_binfhescheme->EvalBinGate(m_params, gate, m_BTKey, ctvector);
}

LWECiphertext BinFHEContext::Bootstrap(ConstLWECiphertext& ct) const {
    return m_binfhescheme->Bootstrap(m_params, m_BTKey, ct);
}

LWECiphertext BinFHEContext::EvalNOT(ConstLWECiphertext& ct) const {
    return m_binfhescheme->EvalNOT(m_params, ct);
}

LWECiphertext BinFHEContext::EvalConstant(bool value) const {
    return m_LWEscheme->NoiselessEmbedding(m_params->GetLWEParams(), value);
}

LWECiphertext BinFHEContext::EvalFunc(ConstLWECiphertext& ct, const std::vector<NativeInteger>& LUT) const {
    return m_binfhescheme->EvalFunc(m_params, m_BTKey, ct, LUT, GetBeta());
}

LWECiphertext BinFHEContext::EvalFloor(ConstLWECiphertext& ct, uint32_t roundbits) const {
    //    auto q = m_params->GetLWEParams()->Getq().ConvertToInt();
    //    if (roundbits != 0) {
    //        NativeInteger newp = this->GetMaxPlaintextSpace();
    //        SetQ(q / newp * (1 << roundbits));
    //    }
    //    SetQ(q);
    //    return res;
    return m_binfhescheme->EvalFloor(m_params, m_BTKey, ct, GetBeta(), roundbits);
}

LWECiphertext BinFHEContext::EvalSign(ConstLWECiphertext& ct, bool schemeSwitch) {
    const auto& params = std::make_shared<BinFHECryptoParams>(*m_params);
    return m_binfhescheme->EvalSign(params, m_BTKey_map, ct, GetBeta(), schemeSwitch);
}

std::vector<LWECiphertext> BinFHEContext::EvalDecomp(ConstLWECiphertext& ct) {
    return m_binfhescheme->EvalDecomp(m_params, m_BTKey_map, ct, GetBeta());
}

std::vector<NativeInteger> BinFHEContext::GenerateLUTviaFunction(NativeInteger (*f)(NativeInteger m, NativeInteger p),
                                                                 NativeInteger p) {
    if (ceil(log2(p.ConvertToInt())) != floor(log2(p.ConvertToInt()))) {
        std::string errMsg("ERROR: Only support plaintext space to be power-of-two.");
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    NativeInteger q        = GetParams()->GetLWEParams()->Getq();
    NativeInteger interval = q / p;
    NativeInteger outerval = interval;
    usint vecSize          = q.ConvertToInt();
    std::vector<NativeInteger> vec(vecSize);
    for (size_t i = 0; i < vecSize; ++i) {
        auto temp = f(NativeInteger(i) / interval, p);
        if (temp >= p) {
            std::string errMsg("ERROR: input function should output in Z_{p_output}.");
            OPENFHE_THROW(not_implemented_error, errMsg);
        }
        vec[i] = temp * outerval;
    }

    return vec;
}




void BinFHEContext::NBTKeyGen(ConstLWEPrivateKey& sk, KEYGEN_MODE keygenMode) {
    //std::cout<<"xzddf btkeygen in binfhecontext.cpp"<<std::endl;
    auto& VNTRUParams = m_params->GetVectorNTRUParams();

    auto temp = VNTRUParams->GetBaseG();

    //预先生成一些密钥
    if (m_timeOptimization) {
        auto gpowermap = VNTRUParams->GetGPowerMap();
        for (std::map<uint32_t, std::vector<NativeInteger>>::iterator it = gpowermap.begin(); it != gpowermap.end();
             ++it) {
            VNTRUParams->Change_BaseG(it->first);
            m_NBTKey_map[it->first] = m_binfhescheme->NKeyGen(m_params, sk, keygenMode);
            // m_NBTKey_map[it->first] = m_binfhescheme->NKeyGen(m_params, sk, keygenMode);
        }
        VNTRUParams->Change_BaseG(temp);
    }

    
    
    //if (m_NBTKey_map.size() != 0) {
    //    m_NBTKey = m_NBTKey_map[temp];
    //}
    //else {
        m_NBTKey           = m_binfhescheme->NKeyGen(m_params, sk, keygenMode);
        // m_NBTKey           = m_binfhescheme->NKeyGen(m_params, sk, keygenMode);
        m_NBTKey_map[temp] = m_NBTKey;
    //}
}


void BinFHEContext::MGNBTKeyGen(std::vector<std::vector<LWEPrivateKey>>& sk, int groupidx, int partyidx, KEYGEN_MODE keygenMode) {
    auto& VNTRUParams = m_params->GetVectorNTRUParams();

    auto temp = VNTRUParams->GetBaseG();

    //预先生成一些密钥
    if (m_timeOptimization) {
        auto gpowermap = VNTRUParams->GetGPowerMap();
        for (std::map<uint32_t, std::vector<NativeInteger>>::iterator it = gpowermap.begin(); it != gpowermap.end();
             ++it) {
            VNTRUParams->Change_BaseG(it->first);
            m_binfhescheme->MGNKeyGen(m_params, sk, groupidx, partyidx, keygenMode);
            // m_NBTKey_map[it->first] = m_binfhescheme->NKeyGen(m_params, sk, keygenMode);
        }
        VNTRUParams->Change_BaseG(temp);
    }

    
    
    //if (m_NBTKey_map.size() != 0) {
    //    m_NBTKey = m_NBTKey_map[temp];
    //}
    //else {
        m_binfhescheme->MGNKeyGen(m_params, sk, groupidx, partyidx, keygenMode);
        // m_NBTKey           = m_binfhescheme->NKeyGen(m_params, sk, keygenMode);
        m_NBTKey_map[temp] = m_NBTKey;
    //}
}


void BinFHEContext::inituser(int groupCount, int groupSize, const std::string& baseDir) {
    namespace fs = std::filesystem;
    fs::path userDir(baseDir);
    if (fs::exists(userDir)) {
        fs::remove_all(userDir);
    }
    fs::create_directories(userDir);
    for (int i = 0; i < groupCount; ++i) {
        fs::path groupPath = userDir / ("group_" + std::to_string(i));
        fs::create_directories(groupPath);
        for (int j = 0; j < groupSize; ++j) {
            fs::path partyPath = groupPath / ("party_" + std::to_string(j));
            fs::create_directories(partyPath);
            // 创建密钥文件
            std::ofstream(partyPath / "sk").close();
            std::ofstream(partyPath / "pk").close();
            std::ofstream(partyPath / "skN").close();
            std::ofstream(partyPath / "invskN").close();
            // 可根据需要添加更多密钥文件
        }
    }
}

void BinFHEContext::LWECRMGen(usint size, const NativeInteger& modulus, int groupIdx) const {
    m_LWEscheme->LWECRMGen(size, modulus, groupIdx);
}

void BinFHEContext::RLWECRPGen(const std::shared_ptr<VectorNTRUCryptoParams>& VNTRUParams) const {
    m_LWEscheme->RLWECRPGen(VNTRUParams);
}



LWEPrivateKey BinFHEContext::LWESecretKeyGen(int groupIdx, int partyIdx) const {
    auto& LWEParams = m_params->GetLWEParams();
    if (LWEParams->GetKeyDist() == GAUSSIAN)
        return m_LWEscheme->LWESecretKeyGenGaussian(LWEParams->Getn(), LWEParams->Getq(), groupIdx, partyIdx);
    return m_LWEscheme->LWESecretKeyGen(LWEParams->Getn(), LWEParams->Getq(), groupIdx, partyIdx);
}

LWEPublicKey BinFHEContext::LWEPublicKeyGen(ConstLWEPrivateKey& sk, int groupIdx, int partyIdx) const {
    return m_LWEscheme->LWEPublicKeyGen(m_params->GetLWEParams(), sk, groupIdx, partyIdx);
}

LWEPublicKey BinFHEContext::LWEJointPublicKeyGen(const std::shared_ptr<LWECryptoParams>& params, int groupIdx) const {
    return m_LWEscheme->LWEJointPublicKeyGen(params, groupIdx);
}

RLWEPrivateKey BinFHEContext::RLWESecretKeyGen(int groupIdx, int partyIdx) const {

    auto& VNTRUParams = m_params->GetVectorNTRUParams();
    NativePoly rlweskPoly(VNTRUParams->GetPolyParams(), Format::COEFFICIENT, true);
    TernaryUniformGeneratorImpl<NativeVector> tug;
    // 生成私钥
    NativeVector rlwesk_vec = tug.GenerateVector(VNTRUParams->GetN(), VNTRUParams->GetQ());
    rlweskPoly.SetValues(rlwesk_vec, Format::COEFFICIENT);
    RLWEPrivateKey rlwesk = std::make_shared<RLWEPrivateKeyImpl>(rlweskPoly);


    std::string rlwesk_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupIdx) + "/party_" + std::to_string(partyIdx);
    try{
        namespace fs = std::filesystem;
        fs::create_directories(rlwesk_dir);
        std::ofstream rlwesk_ofs(rlwesk_dir + "/rlwesk");
        cereal::JSONOutputArchive archive(rlwesk_ofs);
        archive(rlwesk);

    }catch(const std::exception& e){
        std::cerr << "rlwesk save failed. : " << e.what() << std::endl;
    }


    return rlwesk;
}

RLWEPublicKey BinFHEContext::RLWEPublicKeyGen(ConstRLWEPrivateKey& sk, int groupIdx, int partyIdx) const {

    auto& VNTRUParams = m_params->GetVectorNTRUParams();
    
    namespace fs = std::filesystem;
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    std::string rlwesk_dir = user_dir + "/group_" + std::to_string(groupIdx) + "/party_" + std::to_string(partyIdx);
    
    // 1. 从user目录下读取rlweCRP文件并将内容存放到NativePoly变量
    // 如果rlweCRP文件不存在，会输出错误信息并抛出异常
    NativePoly crp(VNTRUParams->GetPolyParams(), Format::COEFFICIENT, true);
    fs::path crp_path = user_dir + "/rlweCRP";
    if (fs::exists(crp_path)) {
        std::ifstream ifs(crp_path);
        NativeVector crp_vec(VNTRUParams->GetN(), VNTRUParams->GetQ());
        // 只读取第一个多项式（前N个值）
        for (size_t j = 0; j < VNTRUParams->GetN(); ++j) {
            string value;
            ifs >> value;
            crp_vec[j] = NativeInteger(value);
        }
        ifs.close();
        crp.SetValues(crp_vec, Format::COEFFICIENT);
    }
    else {
        std::cerr << "错误：rlweCRP文件不存在于路径" << crp_path << std::endl;
        throw std::runtime_error("rlweCRP文件不存在，无法生成RLWE公钥");
    }
    
    // 2. 直接使用传入的sk参数作为rlwesk，不再从文件读取
    NativePoly rlweskPoly = sk->GetElement();
    

    // 3. 从高斯分布生成噪声e（标准差3.0）
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    double STD_DEV = 3.0; // 高斯分布的标准差
    dgg.SetStd(STD_DEV);
    NativePoly e(VNTRUParams->GetPolyParams(), Format::COEFFICIENT, true);
    NativeVector evec = dgg.GenerateVector(VNTRUParams->GetN(), VNTRUParams->GetQ());
    e.SetValues(evec, Format::COEFFICIENT);
    e.SetFormat(Format::EVALUATION);
    
    // 4. 计算公钥r = -rlweCRP*rlwesk + e
    NativePoly rlwecrp = crp; // 先复制crp
    rlweskPoly.SetFormat(Format::EVALUATION);
    rlwecrp.SetFormat(Format::EVALUATION);
    
    auto rPoly = rlweskPoly * rlwecrp;
    rPoly = rPoly.Negate();
    rPoly += e;
    rPoly.SetFormat(Format::COEFFICIENT);
    
    RLWEPublicKey rlwepk = std::make_shared<RLWEPublicKeyImpl>(std::move(crp),std::move(rPoly));

    try{
        namespace fs = std::filesystem;
        fs::create_directories(rlwesk_dir);
        std::ofstream rlwepk_ofs(rlwesk_dir + "/rlwepk");
        cereal::JSONOutputArchive archive(rlwepk_ofs);
        archive(rlwepk);
    }catch(const std::exception& e){
        std::cerr << "rlwepk save failed. : " << e.what() << std::endl;
    }

    
    return rlwepk;
}

RLWEPublicKey BinFHEContext::RLWEJointPublicKeyGen(const std::shared_ptr<VectorNTRUCryptoParams>& params, int groupIdx) const {
    namespace fs = std::filesystem;
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    std::string group_dir = user_dir + "/group_" + std::to_string(groupIdx);
    
    
    // 检查group目录是否存在
    if (!fs::exists(group_dir)) {
        std::cerr << "错误：组目录不存在于路径" << group_dir << std::endl;
        throw std::runtime_error("组目录不存在，无法生成联合公钥");
    }
    
    // 遍历group目录下的所有party目录
    NativePoly jrlwepkr(params->GetPolyParams(), Format::EVALUATION, true);
    NativePoly jrlwepkp(params->GetPolyParams(), Format::COEFFICIENT, true);
    int partyCount = 0;
    for (const auto& entry : fs::directory_iterator(group_dir)) {
        if (entry.is_directory()) {
            std::string dir_name = entry.path().filename().string();
            
            // 检查是否是party目录（格式为"party_数字"）
            if (dir_name.rfind("party_", 0) == 0) {
                std::string rlwepk_path = entry.path().string() + "/rlwepk";
                try{
                    std::ifstream ifs(rlwepk_path);
                    cereal::JSONInputArchive archive(ifs);
                    RLWEPublicKey party_pk = std::make_shared<RLWEPublicKeyImpl>();
                    archive(party_pk);
                    NativePoly partyr = party_pk->Getr();
                    partyr.SetFormat(Format::EVALUATION);
                    jrlwepkr+=partyr;
                    if(partyCount==0){
                        jrlwepkp=party_pk->Getp();
                    }
                    partyCount++;
                }catch(const std::exception& e){
                    std::cerr << "rlwepk load failed. : " << e.what() << std::endl;
                    partyCount++;
                    continue;
                }

            }
        }
    }

    jrlwepkr.SetFormat(Format::COEFFICIENT);
    jrlwepkp.SetFormat(Format::COEFFICIENT);
    //cout<<"jrlwepkr: "<<jrlwepkr<<endl;
    RLWEPublicKey rlwejointPK = std::make_shared<RLWEPublicKeyImpl>(std::move(jrlwepkp),std::move(jrlwepkr));
    
    if (partyCount == 0) {
        std::cerr << "警告：在组" << groupIdx << "中没有找到有效的rlwepk文件" << std::endl;
        throw std::runtime_error("没有找到有效的party公钥，无法生成联合公钥");
    }
    
    // 将联合公钥写入文件
    try{
        fs::path jointrlwepk_path = group_dir + "/jointrlwepk";
        std::ofstream joint_pk_ofs(jointrlwepk_path);
        cereal::JSONOutputArchive archive(joint_pk_ofs);
        archive(rlwejointPK);
    }catch(const std::exception& e){
        std::cerr << "jointrlwepk save failed. : " << e.what() << std::endl;
    }
    

    return rlwejointPK;
}

RLWECiphertext BinFHEContext::AHE_Encrypt(const std::shared_ptr<BinFHECryptoParams>& params, NativePoly m, int index) const{

    auto polyParams =  params->GetVectorNTRUParams()->GetPolyParams();
    uint32_t N = params->GetVectorNTRUParams()->GetN();
    NativeInteger Q = params->GetVectorNTRUParams()->GetQ();
    NativeInteger t = NativeInteger(8) * NativeInteger(N) * NativeInteger(N);

    RLWEPublicKey jpk;
    std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(index);
    std::string rlwejointpk_path = group_dir + "/jointrlwepk";
    try{
        std::ifstream pk_ifs(rlwejointpk_path);
        cereal::JSONInputArchive archive(pk_ifs);
        archive(jpk);
    }catch(const std::exception& e){
        std::cerr << "Error serializing joint RLWE public key to file: " << e.what() << rlwejointpk_path << std::endl;
    }

    //1.选取临时秘密多项式g
    NativePoly gPoly(polyParams, Format::COEFFICIENT, true);
    TernaryUniformGeneratorImpl<NativeVector> tug;
    NativeVector g_vec = tug.GenerateVector(N, Q);
    gPoly.SetValues(g_vec, Format::COEFFICIENT);
    gPoly.SetFormat(Format::EVALUATION);
    

    //2.计算b = r*g + (Q/t) *m + e1;  a = p * g + e2; 
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    dgg.SetStd(3.0);
    NativeVector e1 = dgg.GenerateVector(N, Q);
    NativeVector e2 = dgg.GenerateVector(N, Q);
    NativePoly e1poly(polyParams, Format::COEFFICIENT, true);
    NativePoly e2poly(polyParams, Format::COEFFICIENT, true);
    e1poly.SetValues(e1, Format::COEFFICIENT);
    e1poly.SetFormat(Format::EVALUATION);
    e2poly.SetValues(e2, Format::COEFFICIENT);
    e2poly.SetFormat(Format::EVALUATION);


    auto r = jpk->Getr();
    //cout<<"r: "<<r<<endl;
    auto p = jpk->Getp();
    r.SetFormat(Format::EVALUATION);
    p.SetFormat(Format::EVALUATION);

    // cout<<"m: "<<m<<endl;
    // cout<<"Q/t: "<<(Q/t)<<endl;
    // cout<<"Q/t*m: "<<(Q/t)*m<<endl;
    // auto n1 = (Q/t)*m;
    // n1.SetFormat(Format::COEFFICIENT);
    // cout<<"n1: "<<n1<<endl;
    NativePoly b = r * gPoly + (Q / t) * m + e1poly;
    NativePoly a = p * gPoly + e2poly;

    //是NTT形式
    std::vector<NativePoly> result(2);
    result[0]=b;
    //cout<<"b: "<<b<<endl;
    result[1]=a;
    //cout<<"a: "<<a<<endl;
    return std::make_shared<RLWECiphertextImpl>(std::move(result));
}

NativePoly BinFHEContext::AHE_Decrypt(const std::shared_ptr<BinFHECryptoParams>& params, RLWECiphertext ct, int index, int numofparty) const{

    auto polyParams =  params->GetVectorNTRUParams()->GetPolyParams();
    uint32_t N = params->GetVectorNTRUParams()->GetN();
    NativeInteger Q = params->GetVectorNTRUParams()->GetQ();
    NativeInteger t = NativeInteger(8) * NativeInteger(N) * NativeInteger(N);

    NativePoly b = ct->GetElements()[0];
    //cout<<"b: "<<b<<endl;
    NativePoly a = ct->GetElements()[1];
    //cout<<"a: "<<a<<endl;

    //读取第index组下所有party的rlwesk并存放进RLWEPrivateKeyVec中
    std::vector<RLWEPrivateKey> RLWEPrivateKeyVec(numofparty, RLWEPrivateKey());
    std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(index);
    
    // 遍历所有party的私钥
    for (int i = 0; i < numofparty; ++i) {
        std::string sk_path = group_dir + "/party_" + std::to_string(i) + "/rlwesk";
        try {
            std::ifstream sk_ifs(sk_path);
            if (!sk_ifs) {
                std::cerr << "Failed to open private key file: " << sk_path << std::endl;
                continue;
            }
            cereal::JSONInputArchive archive(sk_ifs);
            RLWEPrivateKey sk;
            archive(sk);
            RLWEPrivateKeyVec[i] = sk;
        } catch (const std::exception& e) {
            std::cerr << "Error loading RLWE private key for party " << i << ": " << e.what() << std::endl;
        }
    }

    for(int i=0;i<numofparty;i++){
        NativePoly z_i = RLWEPrivateKeyVec[i]->GetElement();
        //cout<<"z_i: "<<z_i<<endl;
        z_i.SetFormat(Format::EVALUATION);
        NativePoly tmp = a * z_i;
        b += tmp;
    }
    b.SetFormat(Format::COEFFICIENT);
    // //打印b的前5个系数
    // cout<<"b: "<<b.GetFormat()<< "[";
    // //testPolyDec.SetFormat(Format::COEFFICIENT);
    // auto testDec_values = b.GetValues();
    // for(int i=0; i<std::min(5, static_cast<int>(testDec_values.GetLength())); i++){
    //     cout<<testDec_values[i];
    //     if(i < std::min(5, static_cast<int>(testDec_values.GetLength())) - 1) cout<<" ";
    // }
    // cout<<"]"<<endl;
    
    // NativePoly m = b;
    // for(int i=0;i<N;i++){
    //     m[i] = m[i].MultiplyAndRound(t, Q);
    // }
    //cout<<"b: "<<b<<endl;
    NativePoly m = b.MultiplyAndRound(t, Q);

    //计算的缺陷，在MAR的实现里面，如果值是负数，那么先取其绝对值取整，然后用Q减去取整的值，在这里会使得应该为0的数变成Q
    for(int i=0;i<N;i++){
        if(m[i] == Q){
            m[i] = 0;
        }
    }
   

    

    return m;

}

RLWECiphertext BinFHEContext::AHE_Add(RLWECiphertext ct1, RLWECiphertext ct2) const{
    std::vector<NativePoly> elements;
    elements.reserve(2);
    elements.push_back(ct1->GetElements()[0] + ct2->GetElements()[0]);
    elements.push_back(ct1->GetElements()[1] + ct2->GetElements()[1]);
    return std::make_shared<RLWECiphertextImpl>(elements);
}

RLWECiphertext BinFHEContext::AHE_ScalarMult(RLWECiphertext ct, NativePoly scalar) const{
    if(scalar.GetFormat()!=Format::EVALUATION){
        scalar.SetFormat(Format::EVALUATION);
    }
    std::vector<NativePoly> elements;
    elements.reserve(2);
    elements.push_back(ct->GetElements()[0] * scalar);
    elements.push_back(ct->GetElements()[1] * scalar);
    return std::make_shared<RLWECiphertextImpl>(elements);
}


}  // namespace lbcrypto
