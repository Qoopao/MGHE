#ifndef __RLWE_PUBLICKEY_FWD_H__
#define __RLWE_PUBLICKEY_FWD_H__

#include <memory>

namespace lbcrypto {

class RLWEPublicKeyImpl;

using RLWEPublicKey      = std::shared_ptr<RLWEPublicKeyImpl>;
using ConstRLWEPublicKey = const std::shared_ptr<const RLWEPublicKeyImpl>;

}  // namespace lbcrypto

#endif  // __RLWE_PUBLICKEY_FWD_H__
