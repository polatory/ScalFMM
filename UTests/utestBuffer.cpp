// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#include "FUTester.hpp"

#include "../Src/Containers/FBufferReader.hpp"
#include "../Src/Containers/FBufferWriter.hpp"



/** this class test the buffers container */
class TestBuffer : public FUTester<TestBuffer> {

        // test size
        void TestWriteRead(){
            FBufferWriter writer;

            const int BytesTested = (sizeof(int)+sizeof(char)+sizeof(double)+sizeof(float));
            const int NbTest = 5;
            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                writer << idxWrite << char(idxWrite) << double(idxWrite) << float(idxWrite);
            }

            assert(writer.getSize() == (NbTest*BytesTested));

            FBufferReader reader(writer.getSize());
            assert(reader.getSize() == writer.getSize());

            memcpy(reader.data(), writer.data(), writer.getSize());
            for(int idxRead = 0 ; idxRead < NbTest ; ++idxRead){
                int intval;
                char charval;
                double doubleval;
                float floatval;
                reader >> intval >> charval >> doubleval >> floatval;

                assert(intval == idxRead);
                assert(charval == char(idxRead));
                assert(doubleval == double(idxRead));
                assert(floatval == float(idxRead));

                assert(reader.tell() == (BytesTested * (idxRead+1)));
            }

            assert(reader.tell() == reader.getSize());
            reader.seek(0);
            assert(reader.tell() == 0);

            for(int idxRead = 0 ; idxRead < NbTest ; ++idxRead){
                assert(reader.FBufferReader::getValue<int>() == idxRead);
                assert(reader.FBufferReader::getValue<char>() == char(idxRead));
                assert(reader.FBufferReader::getValue<double>() == double(idxRead));
                assert(reader.FBufferReader::getValue<float>() == float(idxRead));

                assert(reader.tell() == (BytesTested * (idxRead+1)));
            }

            assert(reader.tell() == reader.getSize());
        }


        void TestWriteAt(){
            FBufferWriter writer;

            const int SizeOfInt = int(sizeof(int));
            const int NbTest = 5;
            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                const int position = writer.getSize();

                assert(position == (NbTest * SizeOfInt * idxWrite) + (idxWrite * SizeOfInt));

                writer.FBufferWriter::write<int>(0);

                assert(writer.getSize() == (NbTest * SizeOfInt * idxWrite) + (idxWrite * SizeOfInt) + SizeOfInt);

                for(int count = 1 ; count <= NbTest ; ++count) writer << count;
                writer.writeAt(position, idxWrite);
            }

            FBufferReader reader(writer.getSize());
            memcpy(reader.data(), writer.data(), writer.getSize());

            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                assert(reader.FBufferReader::getValue<int>() == idxWrite);
                for(int count = 1 ; count <= NbTest ; ++count){
                    assert(reader.FBufferReader::getValue<int>() == count);
                }
            }
        }


        // set test
        void SetTests(){
            AddTest(&TestBuffer::TestWriteRead,"Test Write then Read");
            AddTest(&TestBuffer::TestWriteAt,"Test Write at then Read");
        }
};

// You must do this
TestClass(TestBuffer)
