/*
  Copyright 2013 Statoil ASA.

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


#ifndef COMPLETION_HPP_
#define COMPLETION_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Util/Value.hpp>


namespace Opm {

    class DeckKeyword;
    class DeckRecord;

    class Connection {
    public:
        Connection(int i, int j , int k ,
                   int complnum,
                   double depth,
                   WellCompletion::StateEnum state ,
                   const Value<double>& connectionTransmissibilityFactor,
                   const Value<double>& diameter,
                   const Value<double>& skinFactor,
                   const int satTableId,
                   const WellCompletion::DirectionEnum direction = WellCompletion::DirectionEnum::Z);


        bool sameCoordinate(const Connection& other) const;
        bool sameCoordinate(const int i, const int j, const int k) const;

        int getI() const;
        int getJ() const;
        int getK() const;
        double getConnectionTransmissibilityFactor() const;
        const Value<double>& getConnectionTransmissibilityFactorAsValueObject() const;
        double getDiameter() const;
        double getSkinFactor() const;
        int getSatTableId() const;
        void   fixDefaultIJ(int wellHeadI , int wellHeadJ);
        void   shift_complnum( int );
        bool attachedToSegment() const;

        WellCompletion::DirectionEnum getDirection() const;

        bool operator==( const Connection& ) const;
        bool operator!=( const Connection& ) const;

        // related segment number
        // -1 means the completion is not related to segment
        int segment_number = -1;

        double center_depth;
        int complnum;
        WellCompletion::StateEnum state;
        double wellPi = 1.0;

    private:
        int m_i, m_j, m_k;
        Value<double> m_diameter;
        Value<double> m_connectionTransmissibilityFactor;
        Value<double> m_skinFactor;
        int m_satTableId;
        WellCompletion::DirectionEnum m_direction;
        Value<double> getDiameterAsValueObject() const;
        Value<double> getSkinFactorAsValueObject() const;
    };
}



#endif /* COMPLETION_HPP_ */
