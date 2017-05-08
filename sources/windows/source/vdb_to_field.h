
/**
vdb_to_field v0.9
2014-2017 Alexey Abryutin
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
This file may not be redistributed in whole or significant part.
This notice may not be removed or altered from any source distribution.
*/


#pragma warning(disable:4146)
#pragma warning(disable:4308)

#ifndef __vdb_to_field_h__
#define __vdb_to_field_h__


#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

namespace HDK_Sample {

class vdb_to_field : public GAS_SubSolver
{
public:
	GETSET_DATA_FUNCS_F("vel_threshold",   VelThreshold);
	GETSET_DATA_FUNCS_F("srcmult",   SourceMultiply);
	GETSET_DATA_FUNCS_F("time",      Time);
	GETSET_DATA_FUNCS_F("vec_length_threshold", VecLengthThreshold);

    GETSET_DATA_FUNCS_I("nodedepth", NodeDepth);
    GETSET_DATA_FUNCS_I("substeps",  Substeps);
	GETSET_DATA_FUNCS_I("optype",    Operation);
	GETSET_DATA_FUNCS_I("velOpTypeType", vOperation);
	GETSET_DATA_FUNCS_I("normalize", Normalize);

    GETSET_DATA_FUNCS_I("vdbSampler",    VdbSampler);

    GETSET_DATA_FUNCS_S("vdb_name",    VdbName);
    GETSET_DATA_FUNCS_S("vdb_path",    VdbPath);
    GETSET_DATA_FUNCS_S("points_path", PointsPath);

	GETSET_DATA_FUNCS_S("dataname",  DataName);

    GET_DATA_FUNC_S(GAS_NAME_FIELDDEST, FieldDstName);
    GET_DATA_FUNC_S(GAS_NAME_FIELDSOURCE, FieldSrcName);

protected:
    explicit             vdb_to_field(const SIM_DataFactory *factory);
    virtual             ~vdb_to_field();

    bool                 shouldMultiThread(const SIM_RawField *field) const
                         { return field->field()->numTiles() > 1; }

    virtual bool         solveGasSubclass(SIM_Engine &engine,
                                            SIM_Object *obj,
                                            SIM_Time time,
                                            SIM_Time timestep);



private:
    static const SIM_DopDescription     *getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(vdb_to_field, GAS_SubSolver, "vdb_to_field", getDopDescription() );


};

}

#endif