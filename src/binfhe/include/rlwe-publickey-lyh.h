
#ifndef _RLWE_PUBLICKEY_H_
#define _RLWE_PUBLICKEY_H_

#include "lattice/hal/lat-backend.h"
#include "rlwe-publickey-fwd-lyh.h"
#include "math/math-hal.h"
#include "utils/serializable.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lbcrypto {
/**
 * @brief Class that stores the LWE scheme public key; contains a vector
 */
class RLWEPublicKeyImpl : public Serializable {

private:
    NativePoly p;
    NativePoly r;
    
public:
    RLWEPublicKeyImpl() = default;

    explicit RLWEPublicKeyImpl(const NativePoly& p, const NativePoly& r) : p(p), r(r) {}

    RLWEPublicKeyImpl(RLWEPublicKeyImpl&& rhs) noexcept : p(std::move(rhs.p)), r(std::move(rhs.r)) {}   

    RLWEPublicKeyImpl(const RLWEPublicKeyImpl& rhs) : p(rhs.p), r(rhs.r) {}

    RLWEPublicKeyImpl& operator=(const RLWEPublicKeyImpl& rhs) {
        this->p = rhs.p;
        this->r = rhs.r;
        return *this;
    }

    RLWEPublicKeyImpl& operator=(RLWEPublicKeyImpl&& rhs) noexcept {
        this->p = std::move(rhs.p);
        this->r = std::move(rhs.r);
        return *this;
    }

    const NativePoly& Getp() const {
        return p;
    }

    const NativePoly& Getr() const {
        return r;
    }

    void Setp(const NativePoly& p) {
        this->p = p;
    }

    void Setr(const NativePoly& r) {
        this->r = r;
    }


    const NativeInteger& GetModulus() const {
        return r.GetModulus();
    }

    bool operator==(const RLWEPublicKeyImpl& other) const {
        return (p == other.p) && (r == other.r);
    }

    bool operator!=(const RLWEPublicKeyImpl& other) const {
        return !(*this == other);
    }

    template <class Archive>
    void save(Archive& ar, std::uint32_t const version) const {
        ar(::cereal::make_nvp("p", p));
        ar(::cereal::make_nvp("r", r));
    }

    template <class Archive>
    void load(Archive& ar, std::uint32_t const version) {
        if (version > SerializedVersion()) {
            OPENFHE_THROW(deserialize_error, "serialized object version " + std::to_string(version) +
                                                 " is from a later version of the library");
        }

        ar(::cereal::make_nvp("p", p));
        ar(::cereal::make_nvp("r", r));
    }

    std::string SerializedObjectName() const override {
        return "RLWEPublicKey";
    }
    static uint32_t SerializedVersion() {
        return 1;
    }


};

}  // namespace lbcrypto

#endif  // _RLWE_PUBLICKEY_H_
