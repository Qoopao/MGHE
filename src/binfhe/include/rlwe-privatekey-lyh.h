
#ifndef _RLWE_PRIVATEKEY_H_
#define _RLWE_PRIVATEKEY_H_

#include "lattice/hal/lat-backend.h"
#include "lwe-privatekey-fwd.h"
#include "math/math-hal.h"
#include "utils/serializable.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lbcrypto {
/**
 * @brief Class that stores the RLWE scheme secret key; contains a polynomial
 */
class RLWEPrivateKeyImpl : public Serializable {
private:
    NativePoly m_s{};

public:
    RLWEPrivateKeyImpl() = default;

    explicit RLWEPrivateKeyImpl(const NativePoly& s) : m_s(s) {}

    RLWEPrivateKeyImpl(const RLWEPrivateKeyImpl& rhs) : m_s(rhs.m_s) {}

    RLWEPrivateKeyImpl(RLWEPrivateKeyImpl&& rhs) noexcept : m_s(std::move(rhs.m_s)) {}

    RLWEPrivateKeyImpl& operator=(const RLWEPrivateKeyImpl& rhs) {
        this->m_s = rhs.m_s;
        return *this;
    }

    RLWEPrivateKeyImpl& operator=(RLWEPrivateKeyImpl&& rhs) noexcept {
        this->m_s = std::move(rhs.m_s);
        return *this;
    }

    const NativePoly& GetElement() const {
        return m_s;
    }

    void SetElement(const NativePoly& s) {
        m_s = s;
    }

    uint32_t GetLength() const {
        return m_s.GetLength();
    }

    const NativeInteger& GetModulus() const {
        return m_s.GetModulus();
    }

    bool operator==(const RLWEPrivateKeyImpl& other) const {
        return m_s == other.m_s;
    }

    bool operator!=(const RLWEPrivateKeyImpl& other) const {
        return !(*this == other);
    }

    template <class Archive>
    void save(Archive& ar, std::uint32_t const version) const {
        ar(::cereal::make_nvp("s", m_s));
    }

    template <class Archive>
    void load(Archive& ar, std::uint32_t const version) {
        if (version > SerializedVersion()) {
            OPENFHE_THROW(deserialize_error, "serialized object version " + std::to_string(version) +
                                                 " is from a later version of the library");
        }

        ar(::cereal::make_nvp("s", m_s));
    }

    std::string SerializedObjectName() const override {
        return "RLWEPrivateKey";
    }
    static uint32_t SerializedVersion() {
        return 1;
    }


};

}  // namespace lbcrypto

#endif  // _RLWE_PRIVATEKEY_H_
