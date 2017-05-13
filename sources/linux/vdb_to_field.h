/*
* Copyright (c) 2017
*    Side Effects Software Inc.  All rights reserved.
*
* Redistribution and use of Houdini Development Kit samples in source and
* binary forms, with or without modification, are permitted provided that the
* following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
* 2. The name of Side Effects Software may not be used to endorse or
*    promote products derived from this software without specific prior
*    written permission.
*
* THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
* OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
* OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
* NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
* OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
* EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/


#ifndef __VdbToField_h__
#define __VdbToField_h__

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/Dense.h>

#include <GU/GU_PrimVDB.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_Geometry.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_MatrixField.h>

namespace HDK_Sample {

	struct VdbToFieldData
	{
		int susteps;
		int curSubstep;
		int sclrOp;
		int vecOp;
		int maxSrc;

		float simTime;
		float velThreshold;
		float srcMult;

		int srcUsed;
		int srcTotal;

		UT_String vdbName;

		SIM_ScalarField  *dstscalar;
		SIM_VectorField  *dstvector;

		//    data from sop geometry
		std::vector< const GEO_PrimVDB* >*        geoVdbPntrs;
		std::vector< UT_Matrix4D >*               geoMatrix;
		std::vector< float >*                     geoMult;
		std::vector< UT_Vector3F >*               geoV;

		//    data for sourcing
		std::vector< UT_Vector3F >*               vdbV;
		std::vector< float >*                     vdbMult;
		std::vector< openvdb::math::Transform >*  vdbXform;
		std::vector< UT_BoundingBoxF >*           vdbBbox;

	};



	class VdbToField : public GAS_SubSolver
	{
	public:

		GETSET_DATA_FUNCS_F("vel_threshold", VelThreshold);
		GETSET_DATA_FUNCS_F("srcmult", SourceMultiply);
		GETSET_DATA_FUNCS_F("time", Time);
		GETSET_DATA_FUNCS_F("vec_length_threshold", VecLengthThreshold);
		GETSET_DATA_FUNCS_I("nodedepth", NodeDepth);
		GETSET_DATA_FUNCS_I("substeps", Substeps);
		GETSET_DATA_FUNCS_I("maxsources", MaxSources);
		GETSET_DATA_FUNCS_I("optype", Operation);
		GETSET_DATA_FUNCS_I("velOpTypeType", vOperation);
		GETSET_DATA_FUNCS_I("normalize", Normalize);
		GETSET_DATA_FUNCS_I("vdbSampler", VdbSampler);
		GETSET_DATA_FUNCS_S("vdb_name", VdbName);
		GETSET_DATA_FUNCS_S("vdb_path", VdbPath);
		GETSET_DATA_FUNCS_S("dataname", DataName);

		GET_DATA_FUNC_S(GAS_NAME_FIELDDEST, FieldDstName);
		GET_DATA_FUNC_S(GAS_NAME_FIELDSOURCE, FieldSrcName);

	protected:
		explicit             VdbToField(const SIM_DataFactory *factory);
		virtual             ~VdbToField();


		virtual bool         solveGasSubclass(SIM_Engine &engine,
			SIM_Object *obj,
			SIM_Time time,
			SIM_Time timestep);



	private:
		static const SIM_DopDescription     *getDopDescription();

		void getSopSources(
			SIM_Object       *obj,
			const GU_Detail* srcGdp,
			VdbToFieldData&  sd /* solver data */
		);

		template <typename GridType>
		void getVdbSources(
			SIM_Object *obj,
			std::vector< typename GridType::ConstPtr >& vdbGridVec,
			VdbToFieldData&  sd /* solver data */
		);

		DECLARE_STANDARD_GETCASTTOTYPE();
		DECLARE_DATAFACTORY(VdbToField, GAS_SubSolver, "vdb_to_field_wip", getDopDescription());


	};

} // End HDK_Sample namespace

#endif



