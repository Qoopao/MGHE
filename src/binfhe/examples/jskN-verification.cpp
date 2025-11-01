#include "binfhecontext.h"
#include "lattice/hal/lat-backend.h"
#include "math/hal/nativeintbackend.h"
#include <string>
#include <vector>
#include <iostream>
#include <filesystem>

using namespace lbcrypto;
using namespace std;
namespace fs = std::filesystem;

// 从文件读取多项式值
NativePoly ReadPolyFromFile(const std::string& filePath, const std::shared_ptr<VectorNTRUCryptoParams>& polyParams, NativeInteger Q) {
    std::ifstream ifs(filePath);
    if (!ifs.is_open()) {
        std::cerr << "无法打开文件: " << filePath << std::endl;
        return NativePoly(polyParams->GetPolyParams());
    }
    
    NativeVector vec(polyParams->GetPolyParams()->GetRingDimension(), Q);
    for (uint64_t i = 0; i < polyParams->GetPolyParams()->GetRingDimension(); ++i) {
        std::string val_str;
        if (!(ifs >> val_str)) {
            break;
        }
        vec[i] = NativeInteger(val_str);
    }
    ifs.close();
    
    NativePoly poly(polyParams->GetPolyParams());
    poly.SetValues(vec, Format::COEFFICIENT);
    return poly;
}

// 将多项式值写入文件，用于调试
void WritePolyToFile(const NativePoly& poly, const std::string& filePath) {
    std::ofstream ofs(filePath);
    if (!ofs.is_open()) {
        std::cerr << "无法创建文件: " << filePath << std::endl;
        return;
    }
    
    NativePoly tempPoly = poly;
    tempPoly.SetFormat(Format::COEFFICIENT);
    
    for (uint64_t i = 0; i < tempPoly.GetRingDimension(); ++i) {
        ofs << tempPoly[i];
        if (i + 1 != tempPoly.GetRingDimension()) {
            ofs << " ";
        }
    }
    ofs << '\n';
    ofs.close();
}

// 比较两个多项式是否相等
bool ComparePolynomials(const NativePoly& poly1, const NativePoly& poly2, NativeInteger Q) {
    NativePoly temp1 = poly1;
    NativePoly temp2 = poly2;
    
    temp1.SetFormat(Format::COEFFICIENT);
    temp2.SetFormat(Format::COEFFICIENT);
    
    if (temp1.GetRingDimension() != temp2.GetRingDimension()) {
        return false;
    }
    
    for (uint64_t i = 0; i < temp1.GetRingDimension(); ++i) {
        // 确保比较时考虑模Q
        NativeInteger val1 = temp1[i].Mod(Q);
        NativeInteger val2 = temp2[i].Mod(Q);
        if (val1 != val2) {
            std::cout << "多项式在索引" << i << "处不相等: " << val1 << " vs " << val2 << std::endl;
            return false;
        }
    }
    
    return true;
}

int main() {
    try {
        // 初始化BinFHE上下文
        auto cc = BinFHEContext();
        cc.GenerateBinFHEContext(BIGP, XZDDF);
        
        // 获取参数
        auto params = cc.GetParams();
        auto VNTRUParams = params->GetVectorNTRUParams();
        auto polyParams = VNTRUParams->GetPolyParams();
        NativeInteger Q = VNTRUParams->GetQ();
        //uint64_t N = polyParams->GetRingDimension();
        
        // 设置组数和每组人数
        int numGroups = 1;
        int partiesPerGroup = 2;
        
        // 初始化用户目录
        cc.inituser(numGroups, partiesPerGroup);
        
        // 为每个组和用户生成密钥
        std::vector<std::vector<LWEPrivateKey>> sk(numGroups, std::vector<LWEPrivateKey>(partiesPerGroup));
        for (int i = 0; i < numGroups; ++i) {
            for (int j = 0; j < partiesPerGroup; ++j) {
                sk[i][j] = cc.LWESecretKeyGen(i, j);
                cc.LWEPublicKeyGen(sk[i][j], i, j);
                cc.MGNBTKeyGen(sk, i, j);
            }
        }
        
        std::cout << "已生成所有密钥。\n";
        
        // 调用GroupKeyGen生成组密钥
        // std::cout << "正在生成组密钥...\n";
        // cc.GetBinFHEScheme()->GroupNTRUKeyGen(params, numGroups, partiesPerGroup);
        // std::cout << "组密钥生成完成。\n";
        
        // 验证每个组的jskN是否正确
        for (int groupidx = 0; groupidx < numGroups; ++groupidx) {
            std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupidx);
            
            // 检查组目录是否存在
            if (!fs::exists(group_dir) || !fs::is_directory(group_dir)) {
                std::cerr << "组目录不存在: " << group_dir << std::endl;
                continue;
            }
            
            // 读取GroupKeyGen生成的jskN
            std::string jskN_path = group_dir + "/jskN";
            if (!fs::exists(jskN_path) || !fs::is_regular_file(jskN_path)) {
                std::cerr << "jskN文件不存在: " << jskN_path << std::endl;
                continue;
            }
            
            NativePoly jskN_from_file = ReadPolyFromFile(jskN_path, VNTRUParams, Q);
            jskN_from_file.SetFormat(Format::EVALUATION);
                
            // 手动计算所有party的skN的乘积
            NativePoly manual_product(polyParams, Format::COEFFICIENT, true);
            manual_product[0] = NativeInteger(1); // 设置为1
            manual_product.SetFormat(Format::EVALUATION);
            // cout<<"manual_product: "<<endl;
            // for (uint32_t i = 0; i < VNTRUParams->GetN(); ++i) {
            //     cout<<manual_product[i]<<" ";
            // }
            // cout<<endl;
            
            bool hasValidParties = false;
            for (int partyidx = 0; partyidx < partiesPerGroup; ++partyidx) {
                std::string party_dir = group_dir + "/party_" + std::to_string(partyidx);
                std::string skN_path = party_dir + "/skN";
                
                if (fs::exists(skN_path) && fs::is_regular_file(skN_path)) {
                    NativePoly skN_poly = ReadPolyFromFile(skN_path, VNTRUParams, Q);
                    skN_poly.SetFormat(Format::EVALUATION);
                    manual_product = manual_product * skN_poly;
                    hasValidParties = true;
                }
            }
            
            if (hasValidParties) {
                // 转换为COEFFICIENT格式进行比较
                manual_product.SetFormat(Format::COEFFICIENT);
                jskN_from_file.SetFormat(Format::COEFFICIENT);
                
                // 验证结果
                std::cout << "验证组 " << groupidx << " 的jskN..." << std::endl;
                if (ComparePolynomials(manual_product, jskN_from_file, Q)) {
                    std::cout << "✅ 验证成功: jskN等于所有party的skN的乘积。" << std::endl;
                } else {
                    std::cout << "❌ 验证失败: jskN不等于所有party的skN的乘积。" << std::endl;
                    
                    // 保存手动计算的结果用于调试
                    WritePolyToFile(manual_product, group_dir + "/manual_jskN");
                    std::cout << "手动计算的结果已保存到: " << group_dir << "/manual_jskN" << std::endl;
                }
            } else {
                std::cerr << "组 " << groupidx << " 中没有有效的party密钥文件。" << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}