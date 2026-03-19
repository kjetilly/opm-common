/*
  Copyright 2022 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_COPYABLE_PTR_HPP
#define OPM_COPYABLE_PTR_HPP

#include <memory>
#include <opm/common/utility/gpuDecorators.hpp>

namespace Opm {
namespace Utility {

// Wraps a raw pointer and makes it copyable, with GPU support.
//
// On the host: owns the pointed-to object (deep copies on copy-construct/assign,
//   destroys on destruction) — same semantics as the original unique_ptr wrapper.
// On the device (CUDA/HIP kernel): behaves as a non-owning view; copy/assign
//   simply copy the raw pointer and the destructor is a no-op.
//
// WARNING: This template should not be used with polymorphic classes.
//   That would require a virtual clone() method. It will only ever copy
//   the static class type of the pointed-to object.
template <class T>
class CopyablePtr {
public:
    OPM_HOST_DEVICE CopyablePtr() : ptr_(nullptr) {}

    OPM_HOST_DEVICE CopyablePtr(const CopyablePtr& other) {
        if constexpr (OPM_IS_INSIDE_HOST_FUNCTION) {
            ptr_ = other.ptr_ ? new T(*other.ptr_) : nullptr;
        } else {
            ptr_ = other.ptr_; // non-owning on device
        }
    }

    OPM_HOST_DEVICE CopyablePtr(CopyablePtr&& other) noexcept : ptr_(other.ptr_) {
        if constexpr (OPM_IS_INSIDE_HOST_FUNCTION) {
            other.ptr_ = nullptr;
        }
    }

    // copy assignment
    OPM_HOST_DEVICE CopyablePtr& operator=(const CopyablePtr& other) {
        if constexpr (OPM_IS_INSIDE_HOST_FUNCTION) {
            if (this != &other) {
                delete ptr_;
                ptr_ = other.ptr_ ? new T(*other.ptr_) : nullptr;
            }
        } else {
            ptr_ = other.ptr_; // non-owning on device
        }
        return *this;
    }

    // move assignment
    OPM_HOST_DEVICE CopyablePtr& operator=(CopyablePtr&& other) noexcept {
        if constexpr (OPM_IS_INSIDE_HOST_FUNCTION) {
            if (this != &other) {
                delete ptr_;
                ptr_ = other.ptr_;
                other.ptr_ = nullptr;
            }
        } else {
            ptr_ = other.ptr_; // non-owning on device
        }
        return *this;
    }

    // assign directly from a unique_ptr (host only)
    CopyablePtr& operator=(std::unique_ptr<T>&& uptr) {
        delete ptr_;
        ptr_ = uptr.release();
        return *this;
    }

    OPM_HOST_DEVICE ~CopyablePtr() {
        if constexpr (OPM_IS_INSIDE_HOST_FUNCTION) {
            delete ptr_;
        }
        // no-op on device: ownership is not transferred to/from GPU memory
    }

    // member access operator
    OPM_HOST_DEVICE T* operator->() const { return ptr_; }

    // boolean context operator
    OPM_HOST_DEVICE explicit operator bool() const noexcept { return ptr_ != nullptr; }

    // get a raw pointer to the stored value
    OPM_HOST_DEVICE T* get() const { return ptr_; }

    // release ownership (host only — transfers ownership out)
    T* release() {
        T* tmp = ptr_;
        ptr_ = nullptr;
        return tmp;
    }

private:
    T* ptr_;
};

} // namespace Utility
} // namespace Opm
#endif
