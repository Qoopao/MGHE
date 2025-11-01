#ifndef __HPK_FWD_H__
#define __HPK_FWD_H__

#include <memory>

namespace lbcrypto {

class HPKImpl;

using HPK      = std::shared_ptr<HPKImpl>;
using ConstHPK = const std::shared_ptr<const HPKImpl>;

}  // namespace lbcrypto

#endif  // __HPK_FWD_H__
