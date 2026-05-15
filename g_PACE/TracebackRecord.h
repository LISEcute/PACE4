#ifndef TRACEBACK_RECORD_H
#define TRACEBACK_RECORD_H

#include <cstdio>

#include <QDataStream>
#include <QIODevice>
#include <QtGlobal>

// Binary traceback record used by PACE3main.cpp and PACE3statis.cpp.
//
// File format:
//   10 x int16 values
//    5 x float values
//
// Keep this format unchanged to preserve compatibility with old PACE output.
struct TracebackRecord
{
    qint16 ibuf2[10] = {};
    float  buf3[5]  = {};
};

static_assert(sizeof(qint16) == 2, "TracebackRecord expects 16-bit integers");
static_assert(sizeof(float)  == 4, "TracebackRecord expects 32-bit floats");
//--------------------------------------------------------------------------------
inline bool writeTracebackRecord(FILE* file, const TracebackRecord& rec)
{
    if(!file)
        return false;

    const bool okIbuf =
        fwrite(rec.ibuf2, sizeof(qint16), 10, file) == 10;

    const bool okBuf =
        fwrite(rec.buf3, sizeof(float), 5, file) == 5;

    return okIbuf && okBuf;
}
//--------------------------------------------------------------------------------
inline bool readTracebackRecord(QDataStream& in, TracebackRecord& rec)
{
    QIODevice* dev = in.device();
    if(!dev)
        return false;

    constexpr qint64 recordBytes =
        10 * static_cast<qint64>(sizeof(qint16)) +
         5 * static_cast<qint64>(sizeof(float));

    if(dev->bytesAvailable() < recordBytes)
        return false;

    for(int i = 0; i < 10; ++i)
        in >> rec.ibuf2[i];

    for(int i = 0; i < 5; ++i)
        in >> rec.buf3[i];

    return in.status() == QDataStream::Ok;
}

#endif // TRACEBACK_RECORD_H
