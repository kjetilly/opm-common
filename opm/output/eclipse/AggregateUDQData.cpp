/*
  Copyright 2018 Statoil ASA

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

#include <opm/output/eclipse/AggregateUDQData.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/output/eclipse/AggregateGroupData.hpp>
#include <opm/output/eclipse/InteHEAD.hpp>
#include <opm/output/eclipse/UDQDims.hpp>
#include <opm/output/eclipse/VectorItems/intehead.hpp>
#include <opm/output/eclipse/WriteRestartHelpers.hpp>

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQAssign.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQDefine.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQFunctionTable.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQInput.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQParams.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <array>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

// ###########################################################################
// Class Opm::RestartIO::Helpers::AggregateUDQData
// ---------------------------------------------------------------------------

namespace VI = ::Opm::RestartIO::Helpers::VectorItems;

namespace {

    // maximum number of groups
    std::size_t ngmaxz(const std::vector<int>& inteHead)
    {
        return inteHead[20];
    }

    // maximum number of wells
    std::size_t nwmaxz(const std::vector<int>& inteHead)
    {
        return inteHead[163];
    }

    // function to return true if token is a function
    bool isTokenTypeFunc(const Opm::UDQTokenType token)
    {
        return Opm::UDQ::scalarFunc(token)
            || Opm::UDQ::elementalUnaryFunc(token)
            || (token == Opm::UDQTokenType::table_lookup);
    }

    // function to return true if token is a binary operator: type power
    // (exponentiation)
    bool isTokenTypeBinaryPowOp(const Opm::UDQTokenType token)
    {
        return token == Opm::UDQTokenType::binary_op_pow;
    }

    // function to return true if token is a binary operator: type multiply
    // or divide
    bool isTokenTypeBinaryMulDivOp(const Opm::UDQTokenType token)
    {
        bool type = false;

        const auto type_1 = std::array {
            Opm::UDQTokenType::binary_op_div,
            Opm::UDQTokenType::binary_op_mul,
        };

        for (const auto& tok_type : type_1) {
            if (token == tok_type) {
                type = true;
                break;
            }
        }

        return type;
    }

    // function to return true if token is a binary operator: type add or
    // subtract
    bool isTokenTypeBinaryAddSubOp(const Opm::UDQTokenType token)
    {
        bool type = false;

        const auto type_1 = std::array {
            Opm::UDQTokenType::binary_op_add,
            Opm::UDQTokenType::binary_op_sub,
        };

        for (const auto& tok_type : type_1) {
            if (token == tok_type) {
                type = true;
                break;
            }
        }

        return type;
    }

    // function to return true if token is a binary union operator
    bool isTokenTypeBinaryUnionOp(const Opm::UDQTokenType token)
    {
        bool type = false;

        const auto type_1 = std::array {
            Opm::UDQTokenType::binary_op_uadd,
            Opm::UDQTokenType::binary_op_umul,
            Opm::UDQTokenType::binary_op_umin,
            Opm::UDQTokenType::binary_op_umax,
        };

        for (const auto& tok_type : type_1) {
            if (token == tok_type) {
                type = true;
                break;
            }
        }

        return type;
    }

    // function to return true if token is an open or close parenthesis token
    bool isTokenTypeParen(const Opm::UDQTokenType token)
    {
        bool type = false;

        const auto type_1 = std::array {
            Opm::UDQTokenType::open_paren,
            Opm::UDQTokenType::close_paren,
        };

        for (const auto& tok_type : type_1) {
            if (token == tok_type) {
                type = true;
                break;
            }
        }

        return type;
    }

    // A function to return true if the token is an operator
    bool isOperatorToken(const Opm::UDQTokenType token)
    {
        return Opm::UDQ::scalarFunc(token)
            || Opm::UDQ::elementalUnaryFunc(token)
            || Opm::UDQ::binaryFunc(token)
            || Opm::UDQ::setFunc(token);
    }

    // function to return index number of last binary token not inside
    // bracket that is ending the expression
    int numOperators(const std::vector<Opm::UDQToken>& modTokens)
    {
        return std::count_if(modTokens.begin(), modTokens.end(),
                             [](const auto& token)
                             {
                                 return isOperatorToken(token.type())
                                     || isTokenTypeParen(token.type());
                             });
    }

    // function to return the precedence of the current operator/function
    int opFuncPrec(const Opm::UDQTokenType token)
    {
        int prec = 0;
        if (isTokenTypeFunc(token)) prec = 6;
        if (Opm::UDQ::cmpFunc(token)) prec = 5;
        if (isTokenTypeBinaryPowOp(token)) prec = 4;
        if (isTokenTypeBinaryMulDivOp(token)) prec = 3;
        if (isTokenTypeBinaryAddSubOp(token)) prec = 2;
        if (isTokenTypeBinaryUnionOp(token)) prec = 1;
        return prec;
    }

    struct substOuterParentheses
    {
        std::vector<Opm::UDQToken> highestLevOperators;
        std::map<std::size_t, std::vector<Opm::UDQToken>> substitutedTokens;
        int noleadingOpenPar;
        bool leadChangeSign;
    };

    // function to return
    //      a vector of functions and operators at the highest level,
    //      a map of substituted tokens,
    //      the number of leading open_paren that bracket the whole expression,
    //      a logical flag indicating whether there is a leading change of sign in the expression

    substOuterParentheses
    substitute_outer_parenthesis(const std::vector<Opm::UDQToken>& modTokens,
                                 int                               noLeadOpenPar,
                                 bool                              leadChgSgn)
    {
        std::map <std::size_t, std::vector<Opm::UDQToken>> substTok;
        std::vector<Opm::UDQToken> highLevOp;
        std::vector<std::size_t> startParen;
        std::vector<std::size_t> endParen;
        std::size_t level = 0;
        std::size_t search_pos = 0;
        std::size_t subS_max = 0;

        while (search_pos < modTokens.size()) {
            if (modTokens[search_pos].type() == Opm::UDQTokenType::open_paren  && level == 0) {
                startParen.push_back(search_pos);
                ++level;
            }
            else if (modTokens[search_pos].type() == Opm::UDQTokenType::open_paren) {
                ++level;
            }
            else if (modTokens[search_pos].type() == Opm::UDQTokenType::close_paren && level == 1) {
                endParen.push_back(search_pos);
                --level;
            }
            else if (modTokens[search_pos].type() == Opm::UDQTokenType::close_paren) {
                --level;
            }

            ++search_pos;
        }


        //
        //Store all the operators at the highest level
        //include the ecl_expr tokens and replace the content of parentheses with a comp_expr
        if (startParen.size() >= 1) {
            if (startParen[0] > 0) {
                //First store all tokens before the first start_paren
                for (std::size_t i = 0; i < startParen[0]; ++i) {
                    highLevOp.emplace_back(modTokens[i]);
                }
            }

            //
            // Replace content of all parentheses at the highest level by an comp_expr
            // store all tokens including () for all tokens inside a pair of ()
            // also store the tokens between sets of () and at the end of an expression

            for (std::size_t ind = 0; ind < startParen.size();  ++ind) {
                std::vector<Opm::UDQToken> substringToken;
                for (std::size_t i = startParen[ind]; i < endParen[ind]+1; ++i) {
                    substringToken.emplace_back(modTokens[i]);
                }

                // store the content inside the parenthesis
                substTok.emplace(ind, std::move(substringToken));

                //
                // make the vector of high level tokens
                //
                //first add ecl_expr instead of content of (...)

                highLevOp.emplace_back(std::to_string(ind), Opm::UDQTokenType::comp_expr);
                //
                // store all tokens between end_paren before and start_paren after current ()
                subS_max = (ind == startParen.size()-1) ? modTokens.size() : startParen[ind+1];

                if ((endParen[ind] + 1) < subS_max) {
                    for (std::size_t i = endParen[ind] + 1; i < subS_max; ++i) {
                        highLevOp.emplace_back(modTokens[i]);
                    }
                }
            }
        }
        else {
            //
            // treat the case with no ()
            for (std::size_t i = 0; i < modTokens.size(); ++i) {
                highLevOp.emplace_back(modTokens[i]);
            }
        }

        //
        // check if there is a leading minus-sign (change sign)
        if ((modTokens[0].type() == Opm::UDQTokenType::binary_op_sub)) {
            if (startParen.size() > 0) {
                // if followed by start_paren linked to end_paren before end of data
                // set flag and remove from operator list because it is considered as a highest precedence operator
                // unless () go from token 2 two the end of expression
                if ((startParen[0] == 1)  && (endParen[0] < modTokens.size()-1)) {
                    leadChgSgn = true;
                }
            }
            else {
                // set flag and remove from operator list
                leadChgSgn = true;
            }

            if (leadChgSgn) {
                // remove from operator list because it is considered as a highest precedence operator and is
                // therefore not a normal "binary_op_sub" operator
                std::vector<Opm::UDQToken> temp_high_lev_op(highLevOp.begin()+1, highLevOp.end());
                highLevOp = temp_high_lev_op;
            }
        }
        else if (startParen.size() >= 1) {
            //
            // check for leading start_paren combined with end_paren at end of data
            if ((startParen[0] == 0) && (endParen[0] == modTokens.size()-1)) {
                //
                // remove leading and trailing ()
                const std::vector<Opm::UDQToken> modTokens_red(modTokens.begin()+1, modTokens.end()-1);
                noLeadOpenPar += 1;
                //
                // recursive call to itself to re-interpret the token-input

                substOuterParentheses substOpPar =
                    substitute_outer_parenthesis(modTokens_red, noLeadOpenPar, leadChgSgn);

                highLevOp = substOpPar.highestLevOperators;
                substTok = substOpPar.substitutedTokens;
                noLeadOpenPar = substOpPar.noleadingOpenPar;
                leadChgSgn = substOpPar.leadChangeSign;
            }
        }
        // interpretation of token input is completed return resulting object

        return {
            highLevOp,
            substTok,
            noLeadOpenPar,
            leadChgSgn
        };
    }

    // Categorize function in terms of which token-types are used in formula
    //
    // The define_type is (-) the location among a set of tokens of the
    // "top" of the parse tree (AST - abstract syntax tree) i.e. the
    // location of the lowest precedence operator relative to the total set
    // of operators, functions and open-/close - parenthesis
    int define_type(const std::vector<Opm::UDQToken>& tokens)
    {
        int def_type = 0;
        int noLeadOpenPar = 0;
        bool leadChgSgn = false;

        //
        // analyse the expression

        substOuterParentheses expr = substitute_outer_parenthesis(tokens, noLeadOpenPar, leadChgSgn);

        //
        // loop over high level operators to find operator with lowest precedence and highest index

        int curPrec  = 100;
        std::size_t indLowestPrecOper = 0;
        for (std::size_t ind = 0; ind < expr.highestLevOperators.size(); ++ind) {
            if ((expr.highestLevOperators[ind].type() != Opm::UDQTokenType::ecl_expr) &&
                (expr.highestLevOperators[ind].type() != Opm::UDQTokenType::comp_expr) &&
                (expr.highestLevOperators[ind].type() != Opm::UDQTokenType::number))
            {
                const int tmpPrec = opFuncPrec(expr.highestLevOperators[ind].type());
                if (tmpPrec <= curPrec) {
                    curPrec = tmpPrec;
                    indLowestPrecOper = ind;
                }
            }
        }

        //
        // if lowest precedence operator is the first token (and not equal to change sign)
        // NOTE: also for the case with outer () removed
        if ((!expr.leadChangeSign) && (indLowestPrecOper == 0)) {
            // test if operator is a function (precedence = 6)
            if ((curPrec == 6) || (expr.highestLevOperators[indLowestPrecOper].type() == Opm::UDQTokenType::binary_op_sub) ) {
                def_type = -1;
                def_type -= expr.noleadingOpenPar;
            } else {
                // def type is 1 when for all other situations (ecl-expression or number)
                def_type = 1;
            }
        } else {
            //
            // treat cases which start either with (ecl_experessions, open-parenthes or leadChangeSign
            def_type = (expr.leadChangeSign) ?  -1 : 0;
            def_type -= expr.noleadingOpenPar;
            // calculate position of lowest precedence operator
            // account for leading change sign operator
            for (std::size_t ind = 0; ind <= indLowestPrecOper; ++ind) {
                //
                //count operators, including functions and parentheses (not original ecl_experessions)
                if (isOperatorToken(expr.highestLevOperators[ind].type())) {
                    // single operator - subtract one
                    --def_type;
                } else if (expr.highestLevOperators[ind].type() == Opm::UDQTokenType::comp_expr) {
                    // expression in parentheses -  add all operators
                    std::size_t ind_ce = static_cast<std::size_t>(std::stoi(expr.highestLevOperators[ind].str()));
                    auto indSubstTok = expr.substitutedTokens.find(ind_ce);
                    if (indSubstTok != expr.substitutedTokens.end()) {
                        // count the number of operators & parenthes in this sub-expression
                        def_type -= numOperators(indSubstTok->second);
                    } else {
                        const auto msg = fmt::format("Invalid compound expression index {}", ind_ce);
                        Opm::OpmLog::error(msg);
                        throw std::invalid_argument { msg };
                    }
                }
                else if ((expr.highestLevOperators[ind].type() != Opm::UDQTokenType::ecl_expr) &&
                         (expr.highestLevOperators[ind].type() != Opm::UDQTokenType::number))
                {
                    // unknown token - write warning
                    Opm::OpmLog::warning(fmt::format("Unknown tokenType '{}' in define_type()",
                                                     expr.highestLevOperators[ind].str()));
                }
            }
        }

        return def_type;
    }

    std::vector<int>
    ig_phase(const Opm::Schedule&    sched,
             const std::size_t       simStep,
             const std::vector<int>& inteHead)
    {
        std::vector<int> inj_phase(ngmaxz(inteHead), 0);

        auto update_phase = [](const int phase, const int new_phase) {
            if (phase == 0) {
                return new_phase;
            }

            throw std::logic_error {
                "Cannot write restart files with UDA "
                "control on multiple phases in same group"
            };
        };

        for (const auto* group : sched.restart_groups(simStep)) {
            if ((group == nullptr) || !group->isInjectionGroup()) {
                continue;
            }

            auto& int_phase = (group->name() == "FIELD")
                ? inj_phase.back()
                : inj_phase[group->insert_index() - 1];

            int_phase = 0;
            for (const auto& [phase, int_value] : std::array {
                    std::pair {Opm::Phase::OIL,   1},
                    std::pair {Opm::Phase::WATER, 2},
                    std::pair {Opm::Phase::GAS,   3},
                })
            {
                if (! group->hasInjectionControl(phase)) {
                    continue;
                }

                if (group->injectionProperties(phase).uda_phase()) {
                    int_phase = update_phase(int_phase, int_value);
                }
            }
        }

        return inj_phase;
    }

    std::vector<int>
    iuap_data(const Opm::Schedule&                            sched,
              const std::size_t                               simStep,
              const std::vector<Opm::UDQActive::InputRecord>& iuap)
    {
        // Construct the current list of well or group sequence numbers to
        // output the IUAP array.
        std::vector<int> wg_no{};

        for (std::size_t ind = 0; ind < iuap.size(); ++ind) {
            const auto ctrl   = iuap[ind].control;
            const auto wg_key = Opm::UDQ::keyword(ctrl);

            if ((wg_key == Opm::UDAKeyword::WCONPROD) ||
                (wg_key == Opm::UDAKeyword::WCONINJE) ||
                (wg_key == Opm::UDAKeyword::WELTARG))
            {
                wg_no.push_back(sched.getWell(iuap[ind].wgname, simStep).seqIndex());
            }
            else if ((wg_key == Opm::UDAKeyword::GCONPROD) ||
                     (wg_key == Opm::UDAKeyword::GCONINJE))
            {
                if (iuap[ind].wgname != "FIELD") {
                    const auto& group = sched.getGroup(iuap[ind].wgname, simStep);
                    wg_no.push_back(group.insert_index() - 1);
                }
            }
            else {
                const auto msg = fmt::format("Invalid control keyword {} for UDQ {}",
                                             static_cast<int>(ctrl), iuap[ind].udq);

                Opm::OpmLog::error(msg);

                throw std::invalid_argument { msg };
            }
        }

        return wg_no;
    }

    template <typename T>
    std::pair<bool, int>
    findInVector(const std::vector<T>& vecOfElements, const T& element)
    {
        std::pair<bool, int> result;

        // Find given element in vector
        auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

        if (it != vecOfElements.end()) {
            result.second = std::distance(vecOfElements.begin(), it);
            result.first = true;
        }
        else {
            result.first = false;
            result.second = -1;
        }

        return result;
    }

    namespace iUdq {

        Opm::RestartIO::Helpers::WindowedArray<int>
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<int>;
            int nwin = std::max(udqDims[0], 1);
            return WV {
                WV::NumWindows{ static_cast<std::size_t>(nwin) },
                WV::WindowSize{ static_cast<std::size_t>(udqDims[1]) }
            };
        }

        template <class IUDQArray>
        void staticContrib(const Opm::UDQInput& udq_input, IUDQArray& iUdq)
        {
            if (udq_input.is<Opm::UDQDefine>()) {
                const auto& udq_define = udq_input.get<Opm::UDQDefine>();
                const auto& update_status =  udq_define.status();
                const auto& tokens = udq_define.tokens();

                iUdq[0] = (update_status.first == Opm::UDQUpdate::ON)
                    ? 2 : 0;

                iUdq[1] = define_type(tokens);
            }
            else {
                iUdq[0] = iUdq[1] = 0;
            }

            iUdq[2] = udq_input.index.typed_insert_index;
        }

    } // iUdq

    namespace iUad {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<int>>
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<int>;

            auto iuad = std::optional<WV>{};

            if (const auto numIUAD = udqDims[2]; numIUAD > 0) {
                iuad.emplace(WV::NumWindows{ static_cast<std::size_t>(numIUAD) },
                             WV::WindowSize{ static_cast<std::size_t>(udqDims[3]) });
            }

            return iuad;
        }

        template <class IUADArray>
        void staticContrib(const Opm::UDQActive::OutputRecord& udq_record, IUADArray& iUad, int use_cnt_diff)
        {
            iUad[0] = udq_record.uda_code;
            iUad[1] = udq_record.input_index + 1;

            // entry 3  - unknown meaning - value = 1
            iUad[2] = 1;

            iUad[3] = udq_record.use_count;
            iUad[4] = udq_record.use_index + 1 - use_cnt_diff;
        }

    } // iUad

    namespace zUdn {

        Opm::RestartIO::Helpers::WindowedArray<
            Opm::EclIO::PaddedOutputString<8>
        >
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<
                Opm::EclIO::PaddedOutputString<8>>;
            int nwin = std::max(udqDims[0], 1);
            return WV {
                WV::NumWindows{ static_cast<std::size_t>(nwin) },
                WV::WindowSize{ static_cast<std::size_t>(udqDims[4]) }
            };
        }

        template <class zUdnArray>
        void staticContrib(const Opm::UDQInput& udq_input, zUdnArray& zUdn)
        {
            // entry 1 is udq keyword
            zUdn[0] = udq_input.keyword();
            zUdn[1] = udq_input.unit();
        }

    } // zUdn

    namespace zUdl {

        Opm::RestartIO::Helpers::WindowedArray<
            Opm::EclIO::PaddedOutputString<8>
        >
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<
                Opm::EclIO::PaddedOutputString<8>>;

            const int nwin = std::max(udqDims[0], 1);
            return WV {
                WV::NumWindows{ static_cast<std::size_t>(nwin) },
                WV::WindowSize{ static_cast<std::size_t>(udqDims[5]) }
            };
        }

        template <class zUdlArray>
        void staticContrib(const Opm::UDQInput& input, zUdlArray& zUdl)
        {
            // Write out the input formula if key is a DEFINE udq
            if (! input.is<Opm::UDQDefine>()) {
                return;
            }

            const auto l_sstr    = std::string::size_type {8};
            const auto max_l_str = Opm::UDQDims::entriesPerZUDL() * l_sstr;

            const auto& udq_define = input.get<Opm::UDQDefine>();
            const auto& z_data = udq_define.input_string();

            if (z_data.size() > max_l_str) {
                const auto msg =
                    fmt::format(R"(DEFINE expression for UDQ {} is too long.
  Number of characters {} exceeds upper limit of {}.
  Expression: {})",
                                udq_define.keyword(),
                                z_data.size(), max_l_str,
                                z_data);

                throw std::invalid_argument { msg };
            }

            const auto n_sstr = z_data.size() / l_sstr;
            for (auto i = 0*n_sstr; i < n_sstr; ++i) {
                if (i == 0) {
                    auto temp_str = z_data.substr(i * l_sstr, l_sstr);

                    // If first character is a minus sign, change to ~
                    if (temp_str.compare(0, 1, "-") == 0) {
                        temp_str.replace(0, 1, "~");
                    }

                    zUdl[i] = temp_str;
                }
                else {
                    zUdl[i] = z_data.substr(i * l_sstr, l_sstr);
                }
            }

            // Add remainder of last non-zero string
            if ((z_data.size() % l_sstr) > 0) {
                zUdl[n_sstr] = z_data.substr(n_sstr * l_sstr);
            }
        }

    } // zUdl

    namespace iGph {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<int>>
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<int>;

            auto igph = std::optional<WV>{};

            if (const auto numIGPH = udqDims[6]; numIGPH > 0) {
                igph.emplace(WV::NumWindows{ static_cast<std::size_t>(numIGPH) },
                             WV::WindowSize{ static_cast<std::size_t>(1) });
            }

            return igph;
        }

        template <class IGPHArray>
        void staticContrib(const int  inj_phase,
                           IGPHArray& iGph)
        {
            iGph[0] = inj_phase;
        }

    } // iGph

    namespace iUap {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<int>>
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<int>;

            auto iuap = std::optional<WV>{};

            if (const auto numIUAP = udqDims[7]; numIUAP > 0) {
                iuap.emplace(WV::NumWindows{ static_cast<std::size_t>(numIUAP) },
                             WV::WindowSize{ static_cast<std::size_t>(1) });
            }

            return iuap;
        }

        template <class IUAPArray>
        void staticContrib(const int  wg_no,
                           IUAPArray& iUap)
        {
            iUap[0] = wg_no + 1;
        }

    } // iUap

    namespace dUdf {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<double>>
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<double>;

            auto dudf = std::optional<WV>{};

            if (const auto numFieldUDQs = udqDims[12]; numFieldUDQs > 0) {
                dudf.emplace(WV::NumWindows { static_cast<std::size_t>(numFieldUDQs) },
                             WV::WindowSize { static_cast<std::size_t>(1) });
            }

            return dudf;
        }

        template <class DUDFArray>
        void staticContrib(const Opm::UDQState& udq_state,
                           const std::string&   udq,
                           DUDFArray&           dUdf)
        {
            // Set value for group name "FIELD"
            dUdf[0] = udq_state.has(udq)
                ? udq_state.get(udq)
                : Opm::UDQ::restart_default;
        }

    } // dUdf

    namespace dUdg {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<double>>
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<double>;

            auto dudg = std::optional<WV>{};

            if (const auto numGroupUDQs = udqDims[11]; numGroupUDQs > 0) {
                dudg.emplace(WV::NumWindows{ static_cast<std::size_t>(numGroupUDQs) },
                             WV::WindowSize{ static_cast<std::size_t>(udqDims[10]) });
            }

            return dudg;
        }

        template <class DUDGArray>
        void staticContrib(const Opm::UDQState&                  udq_state,
                           const std::vector<const Opm::Group*>& groups,
                           const std::string&                    udq,
                           const std::size_t                     ngmaxz,
                           DUDGArray&                            dUdg)
        {
            for (std::size_t ind = 0; ind < groups.size(); ++ind) {
                const auto* group = groups[ind];

                const auto useDflt = (group == nullptr)
                    || (ind == ngmaxz - 1)
                    || ! udq_state.has_group_var(group->name(), udq);

                dUdg[ind] = useDflt ? Opm::UDQ::restart_default
                    : udq_state.get_group_var(group->name(), udq);
            }
        }

    } // dUdg

    namespace dUdw {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<double>>
        allocate(const std::vector<int>& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<double>;

            auto dudw = std::optional<WV>{};

            if (const auto numWellUDQs = udqDims[9]; numWellUDQs > 0) {
                const auto numWells = std::max(udqDims[8], 1);

                dudw.emplace(WV::NumWindows{ static_cast<std::size_t>(numWellUDQs) },
                             WV::WindowSize{ static_cast<std::size_t>(numWells) });
            }

            return dudw;
        }

        template <class DUDWArray>
        void staticContrib(const Opm::UDQState&            udq_state,
                           const std::vector<std::string>& wells,
                           const std::string               udq,
                           const std::size_t               nwmaxz,
                           DUDWArray&                      dUdw)
        {
            // Initialize array to the default value for the array
            std::fill_n(dUdw.begin(), nwmaxz, Opm::UDQ::restart_default);

            for (std::size_t ind = 0; ind < wells.size(); ++ind) {
                const auto& wname = wells[ind];

                if (udq_state.has_well_var(wname, udq)) {
                    dUdw[ind] = udq_state.get_well_var(wname, udq);
                }
            }
        }
    } // dUdw
}

// ===========================================================================

Opm::RestartIO::Helpers::AggregateUDQData::
AggregateUDQData(const UDQDims& udqDims)
    : iUDQ_ { iUdq::allocate(udqDims.data()) }
    , iUAD_ { iUad::allocate(udqDims.data()) }
    , zUDN_ { zUdn::allocate(udqDims.data()) }
    , zUDL_ { zUdl::allocate(udqDims.data()) }
    , iGPH_ { iGph::allocate(udqDims.data()) }
    , iUAP_ { iUap::allocate(udqDims.data()) }
      // ------------------------------------------------------------
    , dUDF_ { dUdf::allocate(udqDims.data()) }
    , dUDG_ { dUdg::allocate(udqDims.data()) }
    , dUDW_ { dUdw::allocate(udqDims.data()) }
{}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
captureDeclaredUDQData(const Schedule&         sched,
                       const std::size_t       simStep,
                       const UDQState&         udq_state,
                       const std::vector<int>& inteHead)
{
    const auto udqInput = sched.getUDQConfig(simStep).input();

    this->collectUserDefinedQuantities(udqInput, inteHead);

    this->collectUserDefinedArguments(sched, simStep, inteHead);

    if (this->dUDF_.has_value()) {
        this->collectFieldUDQValues(udqInput, udq_state,
                                    inteHead[VI::intehead::NO_FIELD_UDQS]);
    }

    if (this->dUDG_.has_value()) {
        this->collectGroupUDQValues(udqInput, udq_state, ngmaxz(inteHead),
                                    sched.restart_groups(simStep),
                                    inteHead[VI::intehead::NO_GROUP_UDQS]);
    }

    if (this->dUDW_.has_value()) {
        this->collectWellUDQValues(udqInput, udq_state, nwmaxz(inteHead),
                                   sched.wellNames(simStep),
                                   inteHead[VI::intehead::NO_WELL_UDQS]);
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectUserDefinedQuantities(const std::vector<UDQInput>& udqInput,
                             const std::vector<int>&      inteHead)
{
    const auto expectNumUDQ = inteHead[VI::intehead::NO_WELL_UDQS]
        + inteHead[VI::intehead::NO_GROUP_UDQS]
        + inteHead[VI::intehead::NO_FIELD_UDQS];

    int cnt = 0;
    for (const auto& udq_input : udqInput) {
        const auto udq_index = udq_input.index.insert_index;

        auto iudq = this->iUDQ_[udq_index];
        iUdq::staticContrib(udq_input, iudq);

        auto zudn = this->zUDN_[udq_index];
        zUdn::staticContrib(udq_input, zudn);

        auto zudl = this->zUDL_[udq_index];
        zUdl::staticContrib(udq_input, zudl);

        ++cnt;
    }

    if (cnt != expectNumUDQ) {
        OpmLog::error(fmt::format("Inconsistent total number of UDQs: {}, "
                                  "and sum of field, group, "
                                  "and well UDQs: {}", cnt, expectNumUDQ));
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectUserDefinedArguments(const Schedule&         sched,
                            const std::size_t       simStep,
                            const std::vector<int>& inteHead)
{
    const auto& udq_active = sched[simStep].udq_active.get();
    if (! udq_active) {
        return;
    }

    {
        const auto& udq_records = udq_active.iuad();

        int cnt = 0;
        for (std::size_t index = 0; index < udq_records.size(); ++index) {
            const auto& record = udq_records[index];

            const auto wg_key = Opm::UDQ::keyword(record.control);
            if (((wg_key == Opm::UDAKeyword::GCONPROD) ||
                 (wg_key == Opm::UDAKeyword::GCONINJE)) &&
                (record.wg_name() == "FIELD"))
            {
                continue;
            }

            auto iuad = (*this->iUAD_)[cnt];

            const auto use_count_diff = static_cast<int>(index) - cnt;
            iUad::staticContrib(record, iuad, use_count_diff);

            ++cnt;
        }

        if (cnt != inteHead[VI::intehead::NO_IUADS]) {
            OpmLog::error(fmt::format("Inconsistent number of iuad's: {}, "
                                      "number of iuads from intehead {}.",
                                      cnt, inteHead[VI::intehead::NO_IUADS]));
        }
    }

    {
        const auto iuap_vect = iuap_data(sched, simStep, udq_active.iuap());

        if (iuap_vect.size() != static_cast<std::size_t>(inteHead[VI::intehead::NO_IUAPS])) {
            OpmLog::error(fmt::format("Inconsistent number of iuap's: {}, "
                                      "number of iuap's from intehead {}.",
                                      iuap_vect.size(), inteHead[VI::intehead::NO_IUAPS]));
        }

        for (std::size_t index = 0; index < iuap_vect.size(); ++index) {
            auto iuap = (*this->iUAP_)[index];
            iUap::staticContrib(iuap_vect[index], iuap);
        }
    }

    if (inteHead[VI::intehead::NO_IUADS] > 0) {
        const auto phs = ig_phase(sched, simStep, inteHead);

        for (std::size_t index = 0; index < phs.size(); ++index) {
            auto igph = (*this->iGPH_)[index];
            iGph::staticContrib(phs[index], igph);
        }

        if (phs.size() != static_cast<std::size_t>(inteHead[VI::intehead::NGMAXZ])) {
            OpmLog::error(fmt::format("Inconsistent number of igph's: {}, "
                                      "number of igph's from intehead {}",
                                      phs.size(), inteHead[VI::intehead::NGMAXZ]));
        }
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectFieldUDQValues(const std::vector<UDQInput>& udqInput,
                      const UDQState&              udq_state,
                      const int                    expectNumFieldUDQs)
{
    auto ix = std::size_t {0};

    int cnt = 0;
    for (const auto& udq_input : udqInput) {
        if (udq_input.var_type() == UDQVarType::FIELD_VAR) {
            auto dudf = (*this->dUDF_)[ix];

            dUdf::staticContrib(udq_state, udq_input.keyword(), dudf);

            ++ix;
            ++cnt;
        }
    }

    if (cnt != expectNumFieldUDQs) {
        OpmLog::error(fmt::format("Inconsistent number of DUDF elements: {}, "
                                  "expected number of DUDF elements {}.",
                                  cnt, expectNumFieldUDQs));
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectGroupUDQValues(const std::vector<UDQInput>&     udqInput,
                      const UDQState&                  udqState,
                      const std::size_t                ngmax,
                      const std::vector<const Group*>& groups,
                      const int                        expectedNumGroupUDQs)
{
    auto ix = std::size_t{0};

    int cnt = 0;
    for (const auto& udq_input : udqInput) {
        if (udq_input.var_type() == UDQVarType::GROUP_VAR) {
            auto dudg = (*this->dUDG_)[ix];

            dUdg::staticContrib(udqState, groups,
                                udq_input.keyword(),
                                ngmax, dudg);

            ++ix;
            ++cnt;
        }
    }

    if (cnt != expectedNumGroupUDQs) {
        OpmLog::error(fmt::format("Inconsistent number of DUDG elements: {}, "
                                  "expected number of DUDG elements {}.",
                                  cnt, expectedNumGroupUDQs));
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectWellUDQValues(const std::vector<UDQInput>&    udqInput,
                     const UDQState&                 udqState,
                     const std::size_t               nwmax,
                     const std::vector<std::string>& wells,
                     const int                       expectedNumWellUDQs)
{
    auto ix = std::size_t {0};

    int cnt = 0;
    for (const auto& udq_input : udqInput) {
        if (udq_input.var_type() == UDQVarType::WELL_VAR) {
            auto dudw = (*this->dUDW_)[ix];

            dUdw::staticContrib(udqState, wells,
                                udq_input.keyword(),
                                nwmax, dudw);

            ++ix;
            ++cnt;
        }
    }

    if (cnt != expectedNumWellUDQs) {
        OpmLog::error(fmt::format("Inconsistent number of DUDW elements: {}, "
                                  "expected number of DUDW elements {}.",
                                  cnt, expectedNumWellUDQs));
    }
}
