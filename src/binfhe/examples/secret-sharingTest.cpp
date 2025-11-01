#include <bits/stdint-uintn.h>
#include "binfhe-constants.h"
#include "binfhecontext.h"
#include "lattice/hal/lat-backend.h"
#include "lwe-ciphertext-fwd.h"
#include "lwe-privatekey.h"
#include "math/hal/nativeintbackend.h"
#include "utils/inttypes.h"
#include "secret-sharing.h"
#include <cstdint>
#include <iostream>
#include <vector>
#include <string>

using namespace lbcrypto;


// 打印多项式信息的辅助函数
void printPolyInfo(const NativePoly& poly, const std::string& name) {
    std::cout << name << " coefficients: ";
    auto coeffs = poly.GetValues();
    for (size_t i = 0; i < std::min(size_t(5), coeffs.GetLength()); ++i) {
        std::cout << coeffs[i] << " ";
    }
    if (coeffs.GetLength() > 5) {
        std::cout << "...";
    }
    std::cout << std::endl;
}

int main() {

#ifdef NDEBUG
    std::cout << "NDEBUG 已定义 - assert 被禁用" << std::endl;
#else
    std::cout << "NDEBUG 未定义 - assert 正常工作" << std::endl;
#endif

    bool allpass = true;
    std::cout << "开始测试 AdditiveSecretSharing 类..." << std::endl;
    
    // 初始化上下文
    auto cc = BinFHEContext();
    cc.GenerateBinFHEContext(BIGP, XZDDF);
    auto params = cc.GetParams()->GetVectorNTRUParams();
    
    // 参数设置
    int k = 3; // 分片数量
    int K = 1;
    cc.inituser(K, k);

    cc.RLWECRPGen(cc.GetParams()->GetVectorNTRUParams());
    for(int i = 0; i < K; ++i){
        for(int j = 0; j < k; ++j){
            auto rlwesk = cc.RLWESecretKeyGen(i, j);
            auto rlwepk = cc.RLWEPublicKeyGen(rlwesk, i, j);
        }
        auto jointrlwepk = cc.RLWEJointPublicKeyGen(cc.GetParams()->GetVectorNTRUParams(), i); 
    }

    // 创建 AdditiveSecretSharing 实例
    AdditiveSecretSharing secretSharing(params,k);
    
    
    // 生成测试多项式
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(params->GetQ());
    NativePoly secretPoly(dug, params->GetPolyParams(), Format::EVALUATION);

    NativePoly xPoly(dug, params->GetPolyParams(), Format::EVALUATION);
    NativePoly yPoly(dug, params->GetPolyParams(), Format::EVALUATION);
    NativePoly alphaPoly(dug, params->GetPolyParams(), Format::EVALUATION);
    
    std::cout << "\n1. 测试 Split 和 Recover 函数:" << std::endl;
    printPolyInfo(secretPoly, "原始秘密多项式");
    
    // 测试 Split
    auto shares = secretSharing.Split(secretPoly, k);
    std::cout << "生成了 " << shares.size() << " 个分片" << std::endl;
    for (int i = 0; i < k; ++i) {
        printPolyInfo(shares[i], "分片 " + std::to_string(i));
    }
    
    // 测试 Recover
    NativePoly recoveredPoly = secretSharing.Recover(shares);
    printPolyInfo(recoveredPoly, "恢复的多项式");
    
    // 验证恢复是否正确
    if(secretPoly == recoveredPoly){
        std::cout << "✓ Split 和 Recover 测试通过！" << std::endl;
    }else{
        std::cout << "✗ Split 和 Recover 测试失败！" << std::endl;
        allpass = false;  
    }
    
    std::cout << "\n2. 测试 Add 函数:" << std::endl;
    // 生成 x 和 y 的分片
    auto x_shares = secretSharing.Split(xPoly, k);
    auto y_shares = secretSharing.Split(yPoly, k);
    
    // 计算分享的加法
    auto z_shares = secretSharing.Add(x_shares, y_shares);
    
    // 恢复结果并验证
    NativePoly expected_z = xPoly + yPoly;
    NativePoly actual_z = secretSharing.Recover(z_shares);
    
    printPolyInfo(expected_z, "期望的 x+y");
    printPolyInfo(actual_z, "实际计算的 x+y");
    
    if(expected_z == actual_z){
        std::cout << "✓ Add 测试通过！" << std::endl;
    }else{
        std::cout << "✗ Add 测试失败！" << std::endl;
        allpass = false;  
    }
    
    std::cout << "\n3. 测试 ScalarMult 函数:" << std::endl;


    // 计算标量乘法
    auto scalar_z_shares = secretSharing.ScalarMult(alphaPoly, x_shares);
    
    // 恢复结果并验证
    NativePoly expected_scalar_z = alphaPoly * xPoly;
    NativePoly actual_scalar_z = secretSharing.Recover(scalar_z_shares);
    
    printPolyInfo(expected_scalar_z, "期望的 alpha*x");
    printPolyInfo(actual_scalar_z, "实际计算的 alpha*x");
    
    if(expected_scalar_z==actual_scalar_z){
        std::cout << "✓ ScalarMult 测试通过！" << std::endl;
    }else{
        std::cout << "✗ ScalarMult 测试失败！" << std::endl;
        allpass = false;  
    }
    
    std::cout << "\n4. 测试 Mult 函数（使用三元组）:" << std::endl;
    // // 生成测试用的三元组 (a, b, c) 其中 c = a*b
    // NativePoly aPoly(dug, params->GetPolyParams(), Format::EVALUATION);
    // NativePoly bPoly(dug, params->GetPolyParams(), Format::EVALUATION);
    // NativePoly cPoly = aPoly * bPoly;

    
    // // 生成三元组的分片
    // auto a_shares = secretSharing.Split(aPoly, k);
    // auto b_shares = secretSharing.Split(bPoly, k);
    // auto c_shares = secretSharing.Split(cPoly, k);
    secretSharing.TriplesGen(cc, 0, k);
    auto Triples = secretSharing.getTriples();
    auto a_shares = Triples[0][0];
    auto b_shares = Triples[0][1];
    auto c_shares = Triples[0][2];
    
    // 计算乘法
    auto mult_z_shares = secretSharing.Mult(x_shares, y_shares, a_shares, b_shares, c_shares);
    
    // 恢复结果并验证
    NativePoly expected_mult_z = xPoly * yPoly;
    NativePoly actual_mult_z = secretSharing.Recover(mult_z_shares);
    
    printPolyInfo(expected_mult_z, "期望的 x*y");
    printPolyInfo(actual_mult_z, "实际计算的 x*y");
    
    if(expected_mult_z==actual_mult_z){
        std::cout << "✓ Mult 测试通过！" << std::endl;
    }else{
        std::cout << "✗ Mult 测试失败！" << std::endl;
        allpass = false;  
    }
    
    std::cout << "\n5. 测试 Auto 函数（自同构）:" << std::endl;
    // 计算自同构
    uint32_t sigma = 3; // 自同构参数
    auto auto_z_shares = secretSharing.Auto(x_shares, sigma);
    
    // 恢复结果并验证
    NativePoly expected_auto_z = xPoly.AutomorphismTransform(sigma);
    NativePoly actual_auto_z = secretSharing.Recover(auto_z_shares);
    printPolyInfo(expected_auto_z, "期望的自同构结果");
    printPolyInfo(actual_auto_z, "实际计算的自同构结果");
    
    if(expected_auto_z == actual_auto_z){
        std::cout << "✓ Auto 测试通过！" << std::endl;
    }else{
        std::cout << "✗ Auto 测试失败！" << std::endl;
        allpass = false;  
    }
    

    if(allpass){
        std::cout << "\n🎉 所有测试都通过了！" << std::endl;
    }else{
        std::cout << "\n✗ 有测试失败！" << std::endl;
    }

    return 0;
}