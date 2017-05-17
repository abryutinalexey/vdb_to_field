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





#include "vdb_to_field.h"

#include <SOP/SOP_Node.h>

#include <GEO/GEO_PrimPacked.h>
#include <GU/GU_PrimPacked.h>
#include <GA/GA_PrimitiveJSON.h>

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>


#include <GAS/GAS_SubSolver.h>
#include <UT/UT_VoxelArray.h>
#include <OP/OP_Director.h>

#include <SIM/SIM_Engine.h>
#include <SIM/SIM_Options.h>

#include <SIM/SIM_GeometryCopy.h>

using namespace HDK_Sample;




template <typename SamplerType> class CalculateVdbToField
{
public:
	THREADED_METHOD2(
		CalculateVdbToField,
		true,
		calculateScalar,
		SIM_RawField*, dst,
		std::vector< openvdb::FloatGrid::ConstPtr >&, vdb_grid
	)

		THREADED_METHOD2(
			CalculateVdbToField,
			true,
			calculateVector,
			SIM_VectorField*, dst,
			std::vector< openvdb::VectorGrid::ConstPtr >&, vdb_grid
		)

		int   myLength;
	VdbToFieldData sd;// solver data

	bool skipTile(
		UT_VoxelArrayIteratorF& vit,
		std::vector<uint>& vdbPrimNum,
		UT_Vector3F& voxel_size,
		SIM_RawField *dst
	)

	{

		UT_Vector3F voxelPos;
		UT_VoxelTile< fpreal32 >* mtile;
		UT_Vector3F tileRadVec;

		vdbPrimNum.clear();
		mtile = vit.getTile();

		dst->indexToPos(vit.x(), vit.y(), vit.z(), voxelPos);

		tileRadVec = UT_Vector3F(voxel_size.x()*mtile->xres(), voxel_size.y()*mtile->yres(), voxel_size.z()*mtile->zres());

		tileRadVec += voxelPos;

		UT_BoundingBoxF tileBbox(voxelPos, tileRadVec);

		struct BboxPrimNum
		{
			UT_BoundingBoxF bbox;
			float length;
			uint  primNum;
		};

		struct {
			bool operator()(BboxPrimNum a, BboxPrimNum b) const
			{
				return a.length < b.length;
			}
		} customLess;


		std::vector<BboxPrimNum> vdbBboxSorted;
		vdbBboxSorted.resize(sd.vdbBbox->size());

		for (int i = 0; i< sd.vdbBbox->size(); i++)
		{
			vdbBboxSorted[i].bbox = sd.vdbBbox->at(i);
			vdbBboxSorted[i].primNum = i;
			vdbBboxSorted[i].length = (sd.vdbBbox->at(i).center() - tileBbox.center()).length();
		}

		// sorting by distance between source center and tile center
		std::sort(vdbBboxSorted.begin(), vdbBboxSorted.end(), customLess);


		std::vector<uint> tmpPrimNum;
		for (int i = 0; i<vdbBboxSorted.size(); i++)
		{
			if ((tileBbox.getRadius() + vdbBboxSorted[i].bbox.getRadius()) >= vdbBboxSorted[i].length)
				tmpPrimNum.push_back(vdbBboxSorted[i].primNum);
			else break;
		}

		// if tmpPrimNum size bigger than maxSrc than
		// it push not just first maxSrc closest
		// but take every maxSrc

		// for example with maxSrc = 4 tmpPrimNum.size() = 13

		// tmpPrimNum[ 0 1 2 3 4 5 6 7 8 9 10 11 12] length = 13
		// vdbPrimNum[ 0     3     6     9         ] length = 4


		if (tmpPrimNum.size() > sd.maxSrc)
		{
			int idx;
			for (int i = 0; i<sd.maxSrc; i++)
			{
				idx = i*(tmpPrimNum.size() / sd.maxSrc);
				vdbPrimNum.push_back(tmpPrimNum[idx]);
			}
		}
		else
		{
			for (int i = 0; i<tmpPrimNum.size(); i++)
			{
				vdbPrimNum.push_back(tmpPrimNum[i]);
			}
		}

		return vdbPrimNum.size() == 0;
	}



	void calculateScalarPartial(
		SIM_RawField* dst,
		std::vector< openvdb::FloatGrid::ConstPtr >& vdb_grid,
		const UT_JobInfo &info
	);

	void calculateVectorPartial(
		SIM_VectorField* dst,
		std::vector< openvdb::VectorGrid::ConstPtr >& vdb_grid,
		const UT_JobInfo &info
	);



};

// ---------------------------------------   CALCULATE SCALAR    ---------------------------------------------

template <typename SamplerType>
void CalculateVdbToField < SamplerType >::calculateScalarPartial(
	SIM_RawField *dst,
	std::vector< openvdb::FloatGrid::ConstPtr >& vdb_grid,
	const UT_JobInfo &info
)
{

	UT_VoxelArrayIteratorF vit;
	UT_VoxelTileIteratorF  vit_t;
	vit.setArray(dst->fieldNC());

	UT_Interrupt  *boss = UTgetInterrupt();

	UT_Vector3F voxelP;
	UT_Vector3F point_pos;

	UT_Vector3F voxelSize = dst->getVoxelSize();

	vit.setPartialRange(info.job(), info.numJobs());


	// filling sempler vector using vdbXform
	std::vector< openvdb::tools::GridSampler< openvdb::FloatTree, SamplerType > >  sampler;
	for (int i = 0; i<sd.vdbXform->size(); i++)
	{
		openvdb::tools::GridSampler< openvdb::FloatTree, SamplerType >
			interpolator(vdb_grid[i]->constTree(), sd.vdbXform->at(i));
		sampler.push_back(interpolator);
	}


	float interpolatedValue;
	std::vector<uint> vdbPrimNum;

	// iterating over tiles
	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{

		// searching for vdb which intersecting with tile bbox and pushing it to vdbPrimNum
		// skip tile if not found

		if (skipTile(vit, vdbPrimNum, voxelSize, dst))  continue;

		vit_t.setTile(vit);

		// performing computation using different operations
		switch (sd.sclrOp)
		{
		case 0: //-----------------------   COPY   -------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				interpolatedValue = 0;
				dst->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				const openvdb::Vec3R ijk(voxelP.x(), voxelP.y(), voxelP.z());

				// sampling using primitive numbers found previously
				for (int j = 0; j<vdbPrimNum.size(); j++)
				{
					interpolatedValue += sampler[vdbPrimNum[j]].wsSample(ijk)*sd.vdbMult->at(vdbPrimNum[j]);
				}

				vit_t.setValue(interpolatedValue*sd.srcMult);

			}
			break;

		case 1://-----------------------   ADD   -------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				interpolatedValue = 0;
				dst->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				const openvdb::Vec3R ijk(voxelP.x(), voxelP.y(), voxelP.z());
				for (int j = 0; j<vdbPrimNum.size(); j++)
				{
					interpolatedValue += sampler[vdbPrimNum[j]].wsSample(ijk)*sd.vdbMult->at(vdbPrimNum[j]);
				}

				vit_t.setValue(vit_t.getValue() + interpolatedValue*sd.srcMult);
			}
			break;

		case 2://-----------------------   SUBTRACT   -------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				interpolatedValue = 0;
				dst->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				const openvdb::Vec3R ijk(voxelP.x(), voxelP.y(), voxelP.z());
				for (int j = 0; j<vdbPrimNum.size(); j++)
				{
					interpolatedValue += sampler[vdbPrimNum[j]].wsSample(ijk)*sd.vdbMult->at(vdbPrimNum[j]);
				}

				vit_t.setValue(vit_t.getValue() - interpolatedValue*sd.srcMult);
			}
			break;

		case 3: //-----------------------   MULTIPLY   -------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				interpolatedValue = 0;
				dst->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				const openvdb::Vec3R ijk(voxelP.x(), voxelP.y(), voxelP.z());

				for (int j = 0; j<vdbPrimNum.size(); j++)
				{
					interpolatedValue += sampler[vdbPrimNum[j]].wsSample(ijk)*sd.vdbMult->at(vdbPrimNum[j]);
				}

				vit_t.setValue(vit_t.getValue() * interpolatedValue*sd.srcMult);
			}
			break;

		case 4: //-----------------------   DIVIDE   -------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				interpolatedValue = 0;
				dst->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				const openvdb::Vec3R ijk(voxelP.x(), voxelP.y(), voxelP.z());

				for (int j = 0; j<vdbPrimNum.size(); j++)
				{
					interpolatedValue += sampler[vdbPrimNum[j]].wsSample(ijk)*sd.vdbMult->at(vdbPrimNum[j]);
				}

				vit_t.setValue(vit_t.getValue() / (interpolatedValue*sd.srcMult));
			}
			break;

		case 5: //-----------------------   MAXIMUM   -------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{

				interpolatedValue = -1000000;
				dst->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				const openvdb::Vec3R ijk(voxelP.x(), voxelP.y(), voxelP.z());

				for (int j = 0; j<vdbPrimNum.size(); j++)
				{
					interpolatedValue = std::max(sampler[vdbPrimNum[j]].wsSample(ijk)*sd.srcMult*sd.vdbMult->at(vdbPrimNum[j]), interpolatedValue);
				}

				vit_t.setValue(std::max(vit_t.getValue(), interpolatedValue));
			}
			break;


		case 6: //-----------------------   MINIMUN   -------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				interpolatedValue = 100000;
				dst->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				const openvdb::Vec3R ijk(voxelP.x(), voxelP.y(), voxelP.z());

				for (int j = 0; j<vdbPrimNum.size(); j++)
				{
					interpolatedValue = std::max(sampler[vdbPrimNum[j]].wsSample(ijk)*sd.srcMult*sd.vdbMult->at(vdbPrimNum[j]), interpolatedValue);
				}
				vit_t.setValue(std::min(vit_t.getValue(), interpolatedValue));
			}
			break;

		case 7: //-----------------------   TEST   -------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				vit_t.setValue(vdbPrimNum.size());
			}
			break;

		}

	}


}

// ---------------------------------------   CALCULATE VECTOR    ---------------------------------------------

template <typename SamplerType>
void CalculateVdbToField<SamplerType>::calculateVectorPartial(
	SIM_VectorField* dst,
	std::vector< openvdb::VectorGrid::ConstPtr >& vdb_grid,
	const UT_JobInfo &info
)
{

	UT_VoxelArrayIteratorF vit_x;
	UT_VoxelArrayIteratorF vit_y;
	UT_VoxelArrayIteratorF vit_z;

	UT_VoxelTileIteratorF  vit_t;

	SIM_RawField* velx = dst->getField(0);
	SIM_RawField* vely = dst->getField(1);
	SIM_RawField* velz = dst->getField(2);

	UT_VoxelArrayF* velx_rf = velx->fieldNC();
	UT_VoxelArrayF* vely_rf = vely->fieldNC();
	UT_VoxelArrayF* velz_rf = velz->fieldNC();

	vit_x.setArray(velx->fieldNC());
	vit_y.setArray(velx->fieldNC());
	vit_z.setArray(velx->fieldNC());

	UT_Vector3F voxelP, oldVal, newVal, srcPos;
	UT_Vector3F vxlOfst, vxlOfstY, vxlOfstZ;
	UT_Vector3F voxelSize = velx->getVoxelSize();

	UT_Vector3F velv(0, 0, 0);
	UT_Vector3F tmpVecHou(0, 0, 0);
	UT_Vector3F zeroHou(0, 0, 0);
	UT_Vector3F new_pos(0, 0, 0);

	openvdb::Vec3R val(0, 0, 0);
	openvdb::Vec3R tmpVecVdb(0, 0, 0);
	openvdb::Vec3R zeroVdb(0, 0, 0);

	double prevLen, newLen, maxLen;
	float  maxLen2;

	int y_ix, y_iy, y_iz;
	int z_ix, z_iy, z_iz;

	int xrez, yrez, zrez;
	velx->getVoxelRes(xrez, yrez, zrez);

	if (dst->isFaceSampled())
		vxlOfst = voxelSize;
	else
		vxlOfst = UT_Vector3F(0, 0, 0);

	vxlOfstY = UT_Vector3F(vxlOfst.x() / 2, -vxlOfst.y() / 2, 0);
	vxlOfstZ = UT_Vector3F(vxlOfst.x() / 2, 0, -vxlOfst.z() / 2);

	vit_x.setPartialRange(info.job(), info.numJobs());

	int j = 0;

	std::vector<uint> vdbPrimNum;

	std::vector< openvdb::tools::GridSampler< openvdb::Vec3STree, SamplerType> >  sampler;

	for (int i = 0; i<sd.vdbXform->size(); i++)
	{
		openvdb::tools::GridSampler< openvdb::Vec3STree, SamplerType >
			interpolator(vdb_grid[i]->constTree(), sd.vdbXform->at(i));
		sampler.push_back(interpolator);
	}


	// iterating over tiles in destination field x
	for (vit_x.rewind(); !vit_x.atEnd(); vit_x.advanceTile())
	{
		// searching for vdb which intersecting with tile bbox and pushing it to vdbPrimNum
		// skip tile if not found
		if (skipTile(vit_x, vdbPrimNum, voxelSize, velx)) continue;

		vit_t.setTile(vit_x);


		// performing computation using different operations
		switch (sd.vecOp)
		{
		case 0: //----------------------------------  VDB VALUES  ----------------------------------
			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{

				//  searching for destination field y and z voxel indexes

				//                                  +++++++++
				//                                  +       +
				//                    voxel_x  -->  +   z   +
				//                                  +       +
				//                                  +++++++++
				//                                      ^
				//                                      |
				//                                   voxel_y ( y_ix, y_iy, y_iz )


				if (vit_t.x() == xrez) continue;

				velx->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				srcPos = voxelP + UT_Vector3F(0, vxlOfst.y() / 2, vxlOfst.z() / 2);

				y_ix = z_ix = vit_t.x();
				y_iy = z_iy = vit_t.y();
				y_iz = z_iz = vit_t.z();

				oldVal.assign(vit_t.getValue(), vely_rf->getValue(y_ix, y_iy, y_iz), velz_rf->getValue(z_ix, z_iy, z_iz));

				const openvdb::Vec3R ijk(srcPos.x(), srcPos.y(), srcPos.z());

				val = zeroVdb;
				maxLen = 0.0;

				// sampling using primitive numbers found previously
				for (j = 0; j<vdbPrimNum.size(); j++)
				{
					tmpVecVdb = sampler[vdbPrimNum[j]].wsSample(ijk)*sd.vdbMult->at(vdbPrimNum[j]);
					val += tmpVecVdb;
					maxLen = std::max(tmpVecVdb.length(), maxLen);
				}

				val.normalize();
				val *= maxLen;

				newVal.assign(val[0], val[1], val[2]);

				newVal *= sd.srcMult;

				newLen = newVal.length();
				prevLen = oldVal.length();

				if ((prevLen + newLen) < sd.velThreshold)
				{
					//ADD
					newVal += oldVal;
					vit_t.setValue(newVal.x());
					vely_rf->setValue(y_ix, y_iy, y_iz, newVal.y());
					velz_rf->setValue(z_ix, z_iy, z_iz, newVal.z());

				}
				else
				{
					//MAXIMUM
					newVal += oldVal;
					newVal.normalize();
					newVal *= std::max(newLen, prevLen);
					vit_t.setValue(newVal.x());
					vely_rf->setValue(y_ix, y_iy, y_iz, newVal.y());
					velz_rf->setValue(z_ix, z_iy, z_iz, newVal.z());
				}

			}
			break;

		case 1: //---------------------------------- V ATTRIBUTE  ----------------------------------

			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				if (vit_t.x() == xrez) continue;

				velx->indexToPos(vit_t.x(), vit_t.y(), vit_t.z(), voxelP);
				srcPos = voxelP + UT_Vector3F(0, vxlOfst.y() / 2, vxlOfst.z() / 2);

				y_ix = z_ix = vit_t.x();
				y_iy = z_iy = vit_t.y();
				y_iz = z_iz = vit_t.z();

				oldVal.assign(vit_t.getValue(), vely_rf->getValue(y_ix, y_iy, y_iz), velz_rf->getValue(z_ix, z_iy, z_iz));

				const openvdb::Vec3R ijk(srcPos.x(), srcPos.y(), srcPos.z());

				maxLen = 0.0;
				maxLen2 = 0.0;
				velv = zeroHou;

				for (j = 0; j<vdbPrimNum.size(); j++)
				{

					tmpVecVdb = sampler[vdbPrimNum[j]].wsSample(ijk)*sd.vdbMult->at(vdbPrimNum[j]);
					tmpVecHou = sd.vdbV->at(vdbPrimNum[j])*tmpVecVdb.x();
					velv += tmpVecHou;
					maxLen2 = std::max(tmpVecHou.length(), maxLen2);
				}

				velv.normalize();
				velv *= maxLen2;

				newVal = velv;
				newVal *= sd.srcMult;

				newLen = newVal.length();
				prevLen = oldVal.length();

				if ((prevLen + newLen) < sd.velThreshold)
				{
					//ADD
					newVal += oldVal;
					vit_t.setValue(newVal.x());
					vely_rf->setValue(y_ix, y_iy, y_iz, newVal.y());
					velz_rf->setValue(z_ix, z_iy, z_iz, newVal.z());

				}
				else
				{
					//MAXIMUM
					newVal += oldVal;
					newVal.normalize();
					newVal *= std::max(newLen, prevLen);
					vit_t.setValue(newVal.x());
					vely_rf->setValue(y_ix, y_iy, y_iz, newVal.y());
					velz_rf->setValue(z_ix, z_iy, z_iz, newVal.z());
				}

			}

			break;

		case 2: //----------------------------------  TEST  ----------------------------------


			for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance())
			{
				vit_t.setValue(vdbPrimNum.size());
			}

			break;


		} //----------------------------------  END SWITCH  ---------------------------------


	}

}


void
initializeSIM(void *)
{
	IMPLEMENT_DATAFACTORY(VdbToField);
}


VdbToField::VdbToField(const SIM_DataFactory *factory)
	: BaseClass(factory)
{

}

VdbToField::~VdbToField()
{
}


const SIM_DopDescription *VdbToField::getDopDescription()
{
	static PRM_Name     theDstFieldName(GAS_NAME_FIELDDEST, "Dest Field");
	static PRM_Name     theVdbName("vdb_name", "Vdb Name");
	static PRM_Name     theVdbPath("vdb_path", "Vdb Source Path");

	static PRM_Name opTypeName("optype", "Scalar Operation");
	static PRM_Name opTypeItems[] =
	{
		PRM_Name("copy",     "Copy"),
		PRM_Name("add",      "Add"),
		PRM_Name("subtract", "Subtract"),
		PRM_Name("mult",     "Multiply"),
		PRM_Name("divide",   "Divide"),
		PRM_Name("maximum",  "Maximum"),
		PRM_Name("minimum",  "Miminum"),
		PRM_Name("test",      "Test"),
		PRM_Name(0)
	};

	static PRM_Name velOpTypeName("velOpTypeType", "Velocity From");
	static PRM_Name velOpTypeItems[] =
	{
		PRM_Name("vdb_vel", "VDB"),
		PRM_Name("v",       "v"),
		PRM_Name("test",    "Test"),
		PRM_Name(0)
	};

	static PRM_Name vdbSamplerName("vdbSampler", "VDB Sampler");
	static PRM_Name vdbSamplerItems[] =
	{
		PRM_Name("point",     "Point Sampler"),
		PRM_Name("box",       "Box Sampler"),
		PRM_Name("quadratic", "Quadratic Sampler"),
		PRM_Name(0)
	};

	static PRM_Name     srcMultName("srcmult", "Source Multiply");
	static PRM_Default  srcMultDef(1);
	static PRM_Range    srcMultRange(PRM_RANGE_FREE, 0, PRM_RANGE_UI, 10);

	static PRM_Name     velThresholdName("vel_threshold", "Vel Add/Max Threshold");
	static PRM_Default  velThresholdDef(1);
	static PRM_Range    velThresholdRange(PRM_RANGE_FREE, 0, PRM_RANGE_UI, 100);

	static PRM_Name     substepsName("substeps", "Substeps");
	static PRM_Default  substepsDef(1);
	static PRM_Range    substepsRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 20);

	static PRM_Name     maxsourcesName("maxsources", "Max Sources Per Tile");
	static PRM_Default  maxsourcesDef(50);
	static PRM_Range    maxsourcesRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_UI, 100);

	static PRM_Name     timeName("time", "Time");
	static PRM_Default  timeDef(0, "$T");



	static PRM_ChoiceList opTypeMenu((PRM_ChoiceListType)
		(PRM_CHOICELIST_EXCLUSIVE |
			PRM_CHOICELIST_REPLACE),
		opTypeItems);

	static PRM_ChoiceList velOpTypeMenu((PRM_ChoiceListType)
		(PRM_CHOICELIST_EXCLUSIVE |
			PRM_CHOICELIST_REPLACE),
		velOpTypeItems);


	static PRM_ChoiceList vdbSamplerMenu((PRM_ChoiceListType)
		(PRM_CHOICELIST_EXCLUSIVE |
			PRM_CHOICELIST_REPLACE),
		vdbSamplerItems);


	static PRM_Template          theTemplates[] = {
		PRM_Template(PRM_STRING, 1, &theDstFieldName),
		PRM_Template(PRM_STRING, 1, &theVdbName),
		PRM_Template(PRM_STRING,  PRM_TYPE_DYNAMIC_PATH, 1, &theVdbPath,    0, 0, 0, 0, &PRM_SpareData::sopPath),
		PRM_Template(PRM_ORD,    1, &opTypeName,     PRMoneDefaults, &opTypeMenu),

		PRM_Template(PRM_ORD,    1, &velOpTypeName,  PRMzeroDefaults, &velOpTypeMenu),
		PRM_Template(PRM_ORD,    1, &vdbSamplerName, PRMzeroDefaults, &vdbSamplerMenu),

		PRM_Template(PRM_FLT,    1, &srcMultName,   &srcMultDef, 0,   &srcMultRange),
		PRM_Template(PRM_FLT,    1, &velThresholdName,   &velThresholdDef, 0,   &velThresholdRange),
		PRM_Template(PRM_INT,    1, &substepsName,   &substepsDef,   0,  &substepsRange),
		PRM_Template(PRM_INT,    1, &maxsourcesName, &maxsourcesDef, 0,  &maxsourcesRange),
		PRM_Template(PRM_FLT_J,  1, &timeName,  &timeDef),
		PRM_Template()
	};

	static SIM_DopDescription    theDopDescription(
		true,
		"vdb_to_field",
		"Vdb To Field",
		"$OS",
		classname(),
		theTemplates);

	theDopDescription.setDefaultUniqueDataName(1);

	return &theDopDescription;
}



template <typename GridType>
void VdbToField::getVdbSources(
	SIM_Object *obj,
	std::vector< typename GridType::ConstPtr >& vdbGridVec,
	VdbToFieldData&  sd // solver data
)
{

	UT_Matrix4D mtrx;

	for (int i = 0; i<sd.geoVdbPntrs->size(); i++)
	{
		// skip if total multiply ("multiply" * "density,temperature,etc... float attribute") on primitive is equal to 0
		if (sd.geoMult->at(i) == 0.0) continue;

		const GEO_PrimVDB* vbdPrim = sd.geoVdbPntrs->at(i);
		openvdb::GridBase::ConstPtr vdbBaseGrid = vbdPrim->getConstGridPtr();
		typename GridType::ConstPtr vdbGrid = openvdb::gridConstPtrCast< GridType >(vdbBaseGrid);

		if (vdbGrid.get()->empty()) continue;

		sd.srcUsed++;

		vdbGridVec.push_back(vdbGrid);
		openvdb::math::CoordBBox bb;
		vdbGrid.get()->tree().evalLeafBoundingBox(bb);


		// evaluating vdb grid bounding box
		openvdb::CoordBBox box = vdbGrid->evalActiveVoxelBoundingBox();
		openvdb::Coord minmaxIndex[2] = { box.min(),box.max() };
		openvdb::Vec3d minmaxWorld[2] = { openvdb::Vec3f(FLT_MAX,FLT_MAX,FLT_MAX),openvdb::Vec3f(-FLT_MAX,-FLT_MAX,-FLT_MAX) };

		for (int i = 0; i<8; i++)
		{
			openvdb::Vec3f v(minmaxIndex[i & 1].x(), minmaxIndex[(i & 2) >> 1].y(), minmaxIndex[(i & 4) >> 2].z());
			openvdb::math::Vec3d vworld = vdbGrid->indexToWorld(v);

			for (int k = 0; k<3; k++)
			{
				minmaxWorld[0][k] = std::min(minmaxWorld[0][k], vworld[k]);
				minmaxWorld[1][k] = std::max(minmaxWorld[1][k], vworld[k]);
			}
		}

		float mult;

		mtrx = sd.geoMatrix->at(i);
		openvdb::math::Mat4d mat(mtrx.data());

		// test for affinity and invertibility
		// if you push not affinity or not invertibile matrix you will get crash
		// if not affinity or invertibile initializing only with translates

		if (!openvdb::math::isAffine(mat) || !openvdb::math::isInvertible(mat))
		{
			UT_String msg;
			msg.sprintf("Matrix on point #%0d in not Affine or not Invertible", -1);
			this->addError(obj, SIM_MESSAGE, msg, UT_ERROR_WARNING);
			UT_Vector3D pos;
			mtrx.getTranslates(pos);
			mtrx.identity();
			mtrx.setTranslates(pos);
			mat = openvdb::math::Mat4d(mtrx.data());
		}

		openvdb::math::Transform xform = vdbGrid->transform();
		xform.postMult(mat);
		sd.vdbXform->push_back(xform);

		UT_BoundingBoxF bbox(UT_Vector3F(minmaxWorld[0].x(), minmaxWorld[0].y(), minmaxWorld[0].z()),
			UT_Vector3F(minmaxWorld[1].x(), minmaxWorld[1].y(), minmaxWorld[1].z()));

		bbox.transform(mtrx);

		sd.vdbBbox->push_back(bbox);
		sd.vdbMult->push_back(sd.geoMult->at(i));
		sd.vdbV->push_back(sd.geoV->at(i));
	}

}


void
VdbToField::getSopSources(
	SIM_Object *obj,
	const GU_Detail* srcGdp,
	VdbToFieldData&   sd) // solver data
{

	UT_Matrix4D mtrx;
	UT_Vector3F pt_pos;

	const GEO_Primitive*   srcPrim;

	const GU_Detail*       vdbGdp;
	const GEO_Primitive*   vdbGdpPrim;
	GU_ConstDetailHandle   vdbGdh;
	GEO_PrimVDB*           vbdPrim;

	GA_ROAttributeRef  nameAttrRef;
	GA_ROAttributeRef  multAttrRef;
	GA_ROAttributeRef  maxSubstepsAttrRef;

	GA_ROAttributeRef multAllAttrRef, multByNameAttrRef, vAttrRef;

	// itterate over oll primitives in source gdp
	GA_FOR_ALL_PRIMITIVES(srcGdp, srcPrim)
	{

		bool foundVdb = false;
		std::vector<int> primGridsNum;

		//
		//  substeping
		//
		//  if max substeps is < thes total supsteps then
		//  then it push to sd.geoVdbPntrs not first max substeps
		//  but evenly distributed substeps in whole substep range

		//  12 substeps      [0 1 2 3 4 5 6 7 8 9 10 11]
		//  6  max substeps  [0   2   4   6   8   10   ]
		//

		//


		maxSubstepsAttrRef = srcGdp->findAttribute(GA_ATTRIB_PRIMITIVE, "max_substeps");
		if (maxSubstepsAttrRef)
		{
			GA_ROHandleI maxSubstepsAttr(maxSubstepsAttrRef.getAttribute());
			int maxSubsteps = maxSubstepsAttr.get(srcPrim->getMapOffset());

			if (maxSubsteps>0 && maxSubsteps < sd.susteps)
			{
				int div = sd.susteps / maxSubsteps;
				if (sd.curSubstep % div != 0) continue;
			}
		}

		// if source primitive is packed primitive itereting over all primitives inside packed primitice gdp
		if (srcPrim->getTypeId() == GA_PRIMINTERNALSENTINEL)
		{
			const GU_PrimPacked* packedPrim = UTverify_cast<const GU_PrimPacked*>(srcPrim);
			vdbGdh = packedPrim->getPackedDetail();
			GU_DetailHandleAutoReadLock  vdbGdl(vdbGdh);
			vdbGdp = vdbGdl.getGdp();

			// primitive inside packed gdp must contain primitive name attribute
			// if not print warning and skip
			nameAttrRef = vdbGdp->findAttribute(GA_ATTRIB_PRIMITIVE, "name");
			if (nameAttrRef.isInvalid())
			{
				UT_String msg;
				msg.sprintf("\nPacked primitive #%0d doesn't contain inside primitive attribute \"name\" \n", srcPrim->getMapOffset());
				continue;
			}
			GA_ROHandleS nameAttr(nameAttrRef.getAttribute());

			// itereting over all primitives inside packed primitice gdp
			GA_FOR_ALL_PRIMITIVES(vdbGdp, vdbGdpPrim)
			{
				// if primitive is no vdb or name doesn't match name which
				// specified in microsolver parameter "Vdb Name" skip primitive
				if (vdbGdpPrim->getTypeId() != GA_PRIMVDB) continue;
				if (sd.vdbName != nameAttr.get(vdbGdpPrim->getMapOffset())) continue;

				const GEO_PrimVDB* vbdPrim = UTverify_cast<const GEO_PrimVDB*>(vdbGdpPrim);

				if (sd.dstscalar && vbdPrim->getStorageType() != UT_VDB_FLOAT) continue;
				if (sd.dstvector && vbdPrim->getStorageType() != UT_VDB_VEC3F) continue;

				// otherwise push GEO_PrimVDB* to geoVdbPntrs
				sd.geoVdbPntrs->push_back(vbdPrim);
				packedPrim->getFullTransform4(mtrx);
				sd.geoMatrix->push_back(mtrx);

				sd.srcTotal++;

				multAllAttrRef = srcGdp->findAttribute(GA_ATTRIB_PRIMITIVE, "multiply");
				vAttrRef = srcGdp->findAttribute(GA_ATTRIB_PRIMITIVE, "v");
				multByNameAttrRef = srcGdp->findAttribute(GA_ATTRIB_PRIMITIVE, sd.vdbName);
				// multiplier by name means that if houdini primitive has float attribute
				// which match name of that prtimitive it will use it as additional multiplyer

				float totalMult = 1.0;

				if (multAllAttrRef.isValid())
				{
					GA_ROHandleF multAttr(multAllAttrRef);
					totalMult *= multAttr.get(srcPrim->getMapOffset());
				}

				if (multByNameAttrRef.isValid())
				{
					GA_ROHandleF multAttr(multByNameAttrRef);
					totalMult *= multAttr.get(srcPrim->getMapOffset());
				}

				// if v attribute doesn't exist push zero values
				if (vAttrRef.isValid())
				{
					GA_ROHandleV3 vAttr(vAttrRef);
					sd.geoV->push_back(vAttr.get(srcPrim->getMapOffset()));
				}
				else
				{
					sd.geoV->push_back(UT_Vector3F(0, 0, 0));
				}

				sd.geoMult->push_back(totalMult);
			}

		}

		// if source primitive is not packed primitive but vdb
		// preforming same computation as inside packed
		// but asignig identidy matrix not from "fulltransform"
		else if (srcPrim->getTypeId() == GA_PRIMVDB)
		{
			nameAttrRef = srcGdp->findAttribute(GA_ATTRIB_PRIMITIVE, "name");
			if (nameAttrRef.isInvalid())
			{
				UT_String msg;
				msg.sprintf("name Attribute on VDB doesn't exist", srcPrim->getMapOffset());
				this->addError(obj, SIM_MESSAGE, msg, UT_ERROR_WARNING);
				continue;
			}

			GA_ROHandleS nameAttr(nameAttrRef.getAttribute());
			if (sd.vdbName != nameAttr.get(srcPrim->getMapOffset())) continue;

			const GEO_PrimVDB* vbdPrim = UTverify_cast< const GEO_PrimVDB* >(srcPrim);

			if (sd.dstscalar && vbdPrim->getStorageType() != UT_VDB_FLOAT) continue;
			if (sd.dstvector && vbdPrim->getStorageType() != UT_VDB_VEC3F) continue;

			sd.geoVdbPntrs->push_back(vbdPrim);

			sd.srcTotal++;

			UT_Matrix4D mtrx;
			mtrx.identity();

			sd.geoMatrix->push_back(mtrx);

			multAllAttrRef = srcGdp->findAttribute(GA_ATTRIB_PRIMITIVE, "multiply");
			vAttrRef = srcGdp->findAttribute(GA_ATTRIB_PRIMITIVE, "v");
			multByNameAttrRef = srcGdp->findAttribute(GA_ATTRIB_PRIMITIVE, sd.vdbName);

			float totalMult = 1.0;

			if (multAllAttrRef.isValid())
			{
				GA_ROHandleF multAttr(multAllAttrRef);
				totalMult *= multAttr.get(srcPrim->getMapOffset());
			}

			if (multByNameAttrRef.isValid())
			{
				GA_ROHandleF multAttr(multByNameAttrRef);
				totalMult *= multAttr.get(srcPrim->getMapOffset());
			}

			if (vAttrRef.isValid())
			{
				GA_ROHandleV3 vAttr(vAttrRef);
				sd.geoV->push_back(vAttr.get(srcPrim->getMapOffset()));
			}
			else
			{
				sd.geoV->push_back(UT_Vector3F(0, 0, 0));
			}

			sd.geoMult->push_back(totalMult);
		}

	}
}


bool VdbToField::solveGasSubclass(SIM_Engine &engine,
	SIM_Object *obj,
	SIM_Time time,
	SIM_Time timestep)
{


	float srcMult = getSourceMultiply();
	float gtime;
	float stime = getTime();
	int   susteps = getSubsteps();
	float timeStep = timestep;

	if (srcMult == 0)	return true;

	SIM_ScalarField   *dstscalar;
	SIM_VectorField   *dstvector;
	SIM_ScalarField   *mask_field;

	SIM_DataArray     dst;
	SIM_DataArray     mask;

	UT_String         srcPath;
	UT_String         VdbName;
	UT_String         filedName;
	GU_DetailHandle   srcGdh;
	const GU_Detail*  srcGdp;
	SOP_Node*         srcSop;


	getVdbPath(srcPath);
	getVdbName(VdbName);

	getMatchingData(dst, obj, GAS_NAME_FIELDDEST);
	getMatchingData(mask, obj, "source_mask");

	//    data from sop geometry
	std::vector< const GEO_PrimVDB* >        geoVdbPntrs;
	std::vector< UT_Matrix4D >               geoMatrix;
	std::vector< float >                     geoMult;
	std::vector< UT_Vector3F >               geoV;

	//    data for sourcing
	std::vector< UT_Vector3F >               vdbV;
	std::vector< float >                     vdbMult;
	std::vector< openvdb::math::Transform >  vdbXform;
	std::vector< UT_BoundingBoxF >           vdbBbox;


	VdbToFieldData solverData;

	solverData.vdbName = VdbName;

	solverData.vdbV = &vdbV;
	solverData.vdbMult = &vdbMult;
	solverData.vdbXform = &vdbXform;
	solverData.vdbBbox = &vdbBbox;
	solverData.geoVdbPntrs = &geoVdbPntrs;
	solverData.geoMatrix = &geoMatrix;
	solverData.geoMult = &geoMult;
	solverData.geoV = &geoV;

	for (int s = 0; s<susteps; s++)
	{

		solverData.simTime = getTime();
		solverData.susteps = susteps;
		solverData.curSubstep = s;
		solverData.sclrOp = getOperation();
		solverData.vecOp = getvOperation();
		solverData.velThreshold = getVelThreshold();
		solverData.maxSrc = getMaxSources();
		solverData.srcMult = srcMult;
		solverData.srcTotal = 0;
		solverData.srcUsed = 0;

		gtime = stime - (timeStep / susteps)*s;

		srcSop = getSOPNode(srcPath, 1);
		if (!srcSop)
		{
			addError(obj, SIM_MESSAGE, "Invalid Vdb Path", UT_ERROR_WARNING);
			return true;
		}

		srcSop->forceRecook();

		OP_Context  srcContext(gtime);
		srcGdh = srcSop->getCookedGeoHandle(srcContext);
		GU_DetailHandleAutoReadLock  srcGdl(srcGdh);
		srcGdp = srcGdl.getGdp();

		if (!srcGdp)
		{
			addError(obj, SIM_MESSAGE, "Invalid Vdb Source Path", UT_ERROR_WARNING);
			return true;
		}

		if (!dst.entries())
		{
			addError(obj, SIM_MESSAGE, "Fewer fields specified.", UT_ERROR_WARNING);
			return true;
		}

		dstscalar = SIM_DATA_CAST(dst(0), SIM_ScalarField);
		dstvector = SIM_DATA_CAST(dst(0), SIM_VectorField);

		solverData.dstscalar = dstscalar;
		solverData.dstvector = dstvector;

		// filling geo vectors in solverData
		getSopSources(obj, srcGdp, solverData);

		if (dstscalar)
		{
			SIM_RawField *dstfield = dstscalar->getField();

			// filling vdb vectors in solverData
			std::vector< openvdb::FloatGrid::ConstPtr >  vdb_grid;
			getVdbSources< openvdb::FloatGrid >(obj, vdb_grid, solverData);

			CalculateVdbToField< openvdb::tools::PointSampler >      sourcing_point;
			CalculateVdbToField< openvdb::tools::BoxSampler >        sourcing_box;
			CalculateVdbToField< openvdb::tools::QuadraticSampler >  sourcing_quad;

			sourcing_point.sd = solverData;
			sourcing_box.sd = solverData;
			sourcing_quad.sd = solverData;

			switch (getVdbSampler())
			{
			case 0: sourcing_point.calculateScalar(dstfield, vdb_grid);  break;
			case 1: sourcing_box.calculateScalar(dstfield, vdb_grid);    break;
			case 2: sourcing_quad.calculateScalar(dstfield, vdb_grid);   break;
			}
		}


		if (dstvector)
		{

			std::vector< openvdb::VectorGrid::ConstPtr >  vdb_grid;
			getVdbSources< openvdb::VectorGrid >(obj, vdb_grid, solverData);

			CalculateVdbToField< openvdb::tools::StaggeredPointSampler >     sourcing_point;
			CalculateVdbToField< openvdb::tools::StaggeredBoxSampler >       sourcing_box;
			CalculateVdbToField< openvdb::tools::StaggeredQuadraticSampler > sourcing_quad;

			sourcing_point.sd = solverData;
			sourcing_box.sd = solverData;
			sourcing_quad.sd = solverData;

			switch (getVdbSampler())
			{
			case 0: sourcing_point.calculateVector(dstvector, vdb_grid); break;
			case 1: sourcing_box.calculateVector(dstvector, vdb_grid);   break;
			case 2: sourcing_quad.calculateVector(dstvector, vdb_grid);  break;
			}

		}

		UT_String msg;
		msg.sprintf("\nUsing %0d vdb out of %0d  time:  %.6f\n", solverData.srcUsed, solverData.srcTotal, gtime);

		addError(obj, SIM_MESSAGE, msg, UT_ERROR_MESSAGE);

		if (dstscalar)
			dstscalar->pubHandleModification();
		if (dstvector)
			dstvector->pubHandleModification();


		solverData.vdbV->clear();
		solverData.vdbMult->clear();
		solverData.vdbXform->clear();
		solverData.vdbBbox->clear();
		solverData.geoVdbPntrs->clear();
		solverData.geoMatrix->clear();
		solverData.geoMult->clear();
		solverData.geoV->clear();

	}


	return true;
}



