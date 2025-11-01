//hpk的序列化和反序列化

#ifndef OPENFHE_BINFHE_HPK_CUSTOM_SERIAL_H
#define OPENFHE_BINFHE_HPK_CUSTOM_SERIAL_H

#include "lattice/hal/default/poly.h"
#include "math/hal/nativeintbackend.h"
#include "binfhe-base-scheme.h"
#include "utils/exception.h"
#include "utils/serial.h"
#include <fstream>
#include <vector>
#include <string>
#include <cereal/archives/binary.hpp>

namespace lbcrypto {

/**
 * @brief 自定义HPK序列化函数
 * 使用简化的方式将HPK数据序列化到文件
 * @param filePath 输出文件路径
 * @param hpk HPK数据（三维数组格式：[4][digitsG]）
 * @param digitsG digitsG参数值
 * @return 是否成功
 */
inline bool CustomSerializeHPK(
    const std::string& filePath,
    const std::vector<std::vector<NativePoly>>& hpk,
    uint32_t digitsG
) {
    try {
        std::ofstream file(filePath, std::ios::binary);
        if (!file.is_open()) {
            return false;
        }
        
        // 使用cereal库序列化整个hpk对象
        cereal::BinaryOutputArchive archive(file);
        
        // 写入版本号和digitsG
        uint32_t version = 1;
        archive(version, digitsG);
        
        // 序列化hpk数据
        archive(hpk);
        
        file.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error in CustomSerializeHPK: " << e.what() << std::endl;
        return false;
    }
}

/**
 * @brief 自定义HPK反序列化函数
 * 从文件读取HPK数据
 * @param filePath 输入文件路径
 * @param hpk 输出HPK数据
 * @param params VectorNTRU参数
 * @return 是否成功
 */
inline bool CustomDeserializeHPK(
    const std::string& filePath,
    std::vector<std::vector<NativePoly>>& hpk,
    const std::shared_ptr<VectorNTRUCryptoParams>& params
) {
    try {
        std::ifstream file(filePath, std::ios::binary);
        if (!file.is_open()) {
            OPENFHE_THROW(deserialize_error, "Could not open hpk file: " + filePath);
            return false;
        }
        
        // 使用cereal库反序列化hpk对象
        cereal::BinaryInputArchive archive(file);
        
        // 读取版本号和digitsG
        uint32_t version, digitsG;
        archive(version, digitsG);
        
        // 反序列化hpk数据
        archive(hpk);
        
        // 验证多项式格式是否为EVALUATION
        for (int j = 0; j < 4; ++j) {
            for (uint32_t k = 0; k < digitsG; ++k) {
                if (hpk[j][k].GetFormat() != Format::EVALUATION) {
                    // 确保格式正确
                    hpk[j][k].SetFormat(Format::EVALUATION);
                }
            }
        }
        
        file.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error in CustomDeserializeHPK: " << e.what() << std::endl;
        return false;
    }
}

} // namespace lbcrypto

#endif // OPENFHE_BINFHE_HPK_CUSTOM_SERIAL_H