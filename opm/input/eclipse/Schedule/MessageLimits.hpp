/*
  Copyright 2016 Statoil ASA.

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

#ifndef OPM_MESSAGES_HPP
#define OPM_MESSAGES_HPP

namespace Opm {

    class Deck;
    class DeckKeyword;

    class MessageLimits {
    public:
        MessageLimits();
        explicit MessageLimits(const Deck& deck);

        static MessageLimits serializationTestObject();

        ///Get all the value from MESSAGES keyword.
        long long getMessagePrintLimit() const;
        long long getCommentPrintLimit() const;
        long long getWarningPrintLimit() const;
        long long getProblemPrintLimit() const;
        long long getErrorPrintLimit() const;
        long long getBugPrintLimit() const;
        void setMessagePrintLimit(long long value);
        void setCommentPrintLimit(long long value);
        void setWarningPrintLimit(long long value);
        void setProblemPrintLimit(long long value);
        void setErrorPrintLimit(long long value);
        void setBugPrintLimit(long long value);

        long long getMessageStopLimit() const;
        long long getCommentStopLimit() const;
        long long getWarningStopLimit() const;
        long long getProblemStopLimit() const;
        long long getErrorStopLimit() const;
        long long getBugStopLimit() const;
        void setMessageStopLimit(long long value);
        void setCommentStopLimit(long long value);
        void setWarningStopLimit(long long value);
        void setProblemStopLimit(long long value);
        void setErrorStopLimit(long long value);
        void setBugStopLimit(long long value);

        bool operator==(const MessageLimits& data) const;
        void update(const DeckKeyword& keyword);

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(message_print_limit);
            serializer(comment_print_limit);
            serializer(warning_print_limit);
            serializer(problem_print_limit);
            serializer(error_print_limit);
            serializer(bug_print_limit);
            serializer(message_stop_limit);
            serializer(comment_stop_limit);
            serializer(warning_stop_limit);
            serializer(problem_stop_limit);
            serializer(error_stop_limit);
            serializer(bug_stop_limit);
        }

    private:
        long long message_print_limit;
        long long comment_print_limit;
        long long warning_print_limit;
        long long problem_print_limit;
        long long error_print_limit;
        long long bug_print_limit;
        long long message_stop_limit;
        long long comment_stop_limit;
        long long warning_stop_limit;
        long long problem_stop_limit;
        long long error_stop_limit;
        long long bug_stop_limit;
    };
}

#endif
