#ifndef _HPK_H_
#define _HPK_H_

#include "lattice/hal/lat-backend.h"
#include "math/hal/nativeintbackend.h"
#include "math/math-hal.h"
#include "utils/serializable.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lbcrypto {
/**
 * @brief Class that stores the hpk; contains a vector，注意，里面的多项式全都是NTT形式
 */
class HPKImpl : public Serializable {

private:
    NativeInteger Q{};

    uint32_t digitsG{};

    std::vector<std::vector<NativePoly>> hpk{};
    
    
public:
    HPKImpl() = default;

    explicit HPKImpl(const NativeInteger& Q, const uint32_t digitsG) : Q(Q), digitsG(digitsG) {     
        hpk.resize(4, std::vector<NativePoly>(digitsG));
    }

    HPKImpl(HPKImpl&& rhs) noexcept : Q(std::move(rhs.Q)), digitsG(std::move(rhs.digitsG)), hpk(std::move(rhs.hpk)) {}   

    HPKImpl(const HPKImpl& rhs) : Q(rhs.Q), digitsG(rhs.digitsG), hpk(rhs.hpk) {}

    HPKImpl& operator=(const HPKImpl& rhs) {
        this->Q = rhs.Q;
        this->digitsG = rhs.digitsG;
        this->hpk = rhs.hpk;
        return *this;
    }

    HPKImpl& operator=(HPKImpl&& rhs) noexcept {
        this->Q = std::move(rhs.Q);
        this->digitsG = std::move(rhs.digitsG);
        this->hpk = std::move(rhs.hpk);
        return *this;
    }

    const NativeInteger& GetQ() const {
        return Q;   
    }

    void SetQ(const NativeInteger& Q) {
        this->Q = Q;
    }

    const uint32_t& GetDigitsG() const {
        return digitsG;
    }

    void SetDigitsG(const uint32_t digitsG) {
        this->digitsG = digitsG;
    }

    void SetHPK(const uint32_t i, const uint32_t j, const NativePoly& hpk_ij) {
        this->hpk[i][j] = hpk_ij;
    }

    const NativePoly& GetHPK(const uint32_t i, const uint32_t j) const {
        return hpk[i][j];
    }

    std::vector<NativePoly> GetHPK(const uint32_t i) const {
        return hpk[i];
    }

    const NativeInteger& GetModulus() const {
        return Q;
    }

    bool operator==(const HPKImpl& other) const {
        return (Q == other.Q) && (digitsG == other.digitsG) && (hpk == other.hpk);
    }

    bool operator!=(const HPKImpl& other) const {
        return !(*this == other);
    }

    template <class Archive>
    void save(Archive& ar, std::uint32_t const version) const {
        ar(::cereal::make_nvp("Q", Q));
        ar(::cereal::make_nvp("digitsG", digitsG));
        ar(::cereal::make_nvp("hpk", hpk));
    }

    template <class Archive>
    void load(Archive& ar, std::uint32_t const version) {
        if (version > SerializedVersion()) {
            OPENFHE_THROW(deserialize_error, "serialized object version " + std::to_string(version) +
                                                 " is from a later version of the library");
        }

        ar(::cereal::make_nvp("Q", Q));
        ar(::cereal::make_nvp("digitsG", digitsG));
        ar(::cereal::make_nvp("hpk", hpk));
    }

    std::string SerializedObjectName() const override {
        return "HPK";
    }
    static uint32_t SerializedVersion() {
        return 1;
    }


};

}  // namespace lbcrypto

#endif  // _HPK_H_
