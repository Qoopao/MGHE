//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
 * Custom Modifications:
 * - [This code is the implementation of the algorithm in the paper https://eprint.iacr.org/2023/1564]
 * 
 * This modified section follows the terms of the original BSD 2-Clause License.
 * Other modifications are provided under the terms of the BSD 2-Clause License.
 * See the BSD 2-Clause License text below:
 */


//==================================================================================
// Additional BSD License for Custom Modifications:
//
// Copyright (c) 2023 Binwu Xiang,Kaixing Wang and other contributors
//
// All rights reserved.
//
// Author TPOC: wangkaixing22@mails.ucas.ac.cn
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

#include "binfhe-constants.h"
#include "binfhecontext.h"

using namespace lbcrypto;

int main() {
    // Sample Program: Step 1: Set CryptoContext
    auto cc = BinFHEContext();

    cc.GenerateBinFHEContext(BIGP, XZDDF);

    // Sample Program: Step 2: Key Generation LWE的密钥，高斯的
    auto sk = cc.KeyGen();
    // std::cout<<"q="<<sk->GetModulus()<<endl;
    // for(int i=0;i<805;++i){
    // 	std::cout<<sk->GetElement()[i]<<endl;
    // }

    int m0=1;
    int m1=1;
    LWEPlaintext result;


    // Generate the bootstrapping keys (refresh and switching keys)
    std::cout << "Generating the bootstrapping keys..." << std::endl;
    cc.NBTKeyGen(sk);       //NTRU的密钥，即vec-NTRU加密LWE的s
    std::cout << "Completed the key generation." << std::endl;


    // Sample Program: Step 3: Encryption
    auto ct1 = cc.Encrypt(sk, m0);
    auto ct2 = cc.Encrypt(sk, m1);
    LWECiphertext ctAND1;
    

    {
        {
            //构造(5q/8,0)
            auto n = 16;
            NativeVector zero(n,0);
            uint32_t q = 2048;
            zero.SetModulus(q);
            NativeInteger temp_b= 5*q/8;
            LWECiphertext ct_temp = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(zero), temp_b.Mod(q)));
            // cout<<ct_temp->GetA()<<endl;
            // cout<<ct_temp->GetB()<<endl;
            auto ctprep=ct1;
            // the additive homomorphic operation for XOR/NXOR is different from the other gates we compute
            // 2*(ct1 - ct2) mod 4 for XOR, me map 1,2 -> 1 and 3,0 -> 0
        
            //LWEscheme->EvalAddEq(ctprep, ct2);
        
            ctprep->GetA().ModAddEq(ct2->GetA());
            ctprep->GetB().ModAddFastEq(ct2->GetB(), ct1->GetModulus());
            //LWEscheme->EvalSubEq(ct_temp,ctprep);
            ct_temp->GetA().ModSubEq(ctprep->GetA());
            ct_temp->GetB().ModSubFastEq(ctprep->GetB(), ct_temp->GetModulus());
        
            auto ctNANDT=ct_temp;   //FHEW里面的NAND门
            // printA(ctNANDT->GetA(), "ctNANDT");
            // printB(ctNANDT->GetB(), "ctNANDT");
        
            LWEPlaintext rree;
            cc.Decrypt(sk,ctNANDT,&rree);//round(q/4 * x)
            //cc.MGHE_Decrypt(sk,ctNANDT,&rree);
            std::cout << "Result of NAND computation = " << rree << std::endl;
    }
    }

    // Sample Program: Step 4: Evaluation
    std::cout << "Start  the  gate bootstrapping " << std::endl;

    //Notice: We have only made specific modifications for NAND gates, and will add other gates in the future.
    //这一步门计算包括了bootstrapping，
    ctAND1 = cc.EvalBinGate(NAND, ct1, ct2);

    cc.Decrypt(sk, ctAND1, &result);

    std::cout << "Result of encrypted computation of ( "<<m0<<" NAND "<<m1<<" ) = " << result << std::endl;

    return 0;
}
