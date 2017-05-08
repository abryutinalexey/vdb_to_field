
/**
vdb_to_field v0.9
2014-2017 Alexey Abryutin
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
This file may not be redistributed in whole or significant part.
This notice may not be removed or altered from any source distribution.
*/


#include "vdb_to_field.h"

#include <SOP/SOP_Node.h>

#include <GU/GU_PrimVDB.h>
#include <GEO/GEO_PrimPacked.h>
#include <GU/GU_PrimPacked.h>
#include <GA/GA_PrimitiveJSON.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/Dense.h>

#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_MatrixField.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_Geometry.h>
#include <GAS/GAS_SubSolver.h>
#include <UT/UT_VoxelArray.h>
#include <OP/OP_Director.h>

#include <SIM/SIM_Engine.h>
#include <SIM/SIM_Options.h>

#include <SIM/SIM_GeometryCopy.h>

using namespace HDK_Sample;



template <class T1, class T2> class m_gasvdbtofield
{
public:
  THREADED_METHOD6(
    		   m_gasvdbtofield,
			   true,
			   calculateScalar,
			   SIM_RawField*,   dst,
               std::vector<openvdb::FloatGrid::Ptr >&, vdb_grid,
			   std::vector<UT_BoundingBoxF>&, vdb_bbox,
               std::vector<openvdb::math::Transform >&, vdb_xform,
			   float, src_mult,
			   int, operation
		 )


 THREADED_METHOD7(
    	   m_gasvdbtofield,
		   true,
		   calculateVector,
		   SIM_VectorField*, dst,
		   int, idx,
           std::vector< openvdb::VectorGrid::Ptr >&, vdb_grid,
		   std::vector<UT_BoundingBoxF>&, vdb_bbox,
           std::vector<openvdb::math::Transform >&, vdb_xform,
		   float, src_mult,
		   int, operation
		 )

  int   myLength;
  int   vdb_depth;
  bool  normalize;
  float vec_lenght_threshold;

  std::vector<UT_Vector3F >       vdb_vel_in;
  std::vector< float >            vdb_mult;
  std::vector< openvdb::Mat3R >   vdb_vel_xform;


  bool skipTile(UT_VoxelArrayIteratorF& vit,
				std::vector<UT_BoundingBoxF>& vdb_bbox,
				std::vector<uint>& vdb_prim_num,
				UT_Vector3F& voxel_size,
				SIM_RawField *dst,
				int mode
 	      )

  {
      UT_Vector3F voxel_pos;
      UT_VoxelTile< fpreal32 >* mtile;
      bool in_bbox;
      UT_Vector3F tile_radius_vector;
      int j=0;
      vdb_prim_num.clear();
      in_bbox=true;
      mtile = vit.getTile();
      dst->indexToPos( vit.x(), vit.y(), vit.z(), voxel_pos  );

      tile_radius_vector = UT_Vector3F(  voxel_size.x()*mtile->xres(),  voxel_size.y()*mtile->yres(), voxel_size.z()*mtile->zres()  );

      tile_radius_vector += voxel_pos;

	  if(mode ==0)
	  {
		  for(j=0; j<vdb_bbox.size();j++)
		  {
			UT_BoundingBoxF tile_bbox( voxel_pos,  tile_radius_vector );
			if(tile_bbox.computeIntersection(vdb_bbox[j]))
			{
			  in_bbox = false;
			  vdb_prim_num.push_back(j);
			}
		  }
	  }
	  else
	  {
		  for(j=0; j<vdb_bbox.size();j++)
		  {
			UT_BoundingBoxF tile_bbox( voxel_pos,  tile_radius_vector );
			if(tile_bbox.computeIntersection(vdb_bbox[j]))
			{
			  //in_bbox = true;
			  return false;
			}
		  }
	  }

      return in_bbox;

  }

  void calculateScalarPartial( SIM_RawField *dst,
                               std::vector< openvdb::FloatGrid::Ptr >& vdb_grid,
							   std::vector<UT_BoundingBoxF>& vdb_bbox,
                               std::vector<openvdb::math::Transform >& vdb_xform,
							   float src_mult,
							   int operation,
							   const UT_JobInfo &info );

  void calculateVectorPartial( SIM_VectorField* dst_x,
                                int idx,
                                std::vector< openvdb::VectorGrid::Ptr >& vdb_grid,
								std::vector<UT_BoundingBoxF>& vdb_bbox,
                                std::vector<openvdb::math::Transform >& vdb_xform,
								float src_mult,
								int operation,
								const UT_JobInfo &info );



};
  // ---------------------------------------   CALCULATE SCALAR    ---------------------------------------------

template <class T1, class T2> void m_gasvdbtofield<T1, T2>::calculateScalarPartial( SIM_RawField *dst,
				                                std::vector< openvdb::FloatGrid::Ptr >& vdb_grid,
											    std::vector<UT_BoundingBoxF>& vdb_bbox,
				                                std::vector<openvdb::math::Transform >& vdb_xform,
											    float src_mult,
											    int operation,
											    const UT_JobInfo &info )
  {

    UT_VoxelArrayIteratorF vit;
    UT_VoxelTileIteratorF  vit_t;
    vit.setArray( dst->fieldNC());

    UT_Interrupt                *boss = UTgetInterrupt();

    UT_Vector3F voxel_pos;
    UT_Vector3F point_pos;

    int i =0;
    GA_Offset ptoff;

    UT_Vector3F voxel_size = dst->getVoxelSize();

    float tile_radius;

	float vdbm;

    UT_Vector3F pos_pt_f;

    UT_ValArray<GA_Offset> ptoff_array;

    vit.setPartialRange(info.job(), info.numJobs());

    int j=0;


    //std::vector< openvdb::tools::GridSampler< openvdb::FloatTree, openvdb::tools::PointSampler> >  sampler_vec;

	std::vector< openvdb::tools::GridSampler<T1, T2> >  sampler_vec;

    for(i=0; i<vdb_xform.size(); i++)
    {
        openvdb::tools::GridSampler<T1, T2>
        interpolator( vdb_grid[ i]->constTree(), vdb_xform[i] );
        sampler_vec.push_back( interpolator );
    }


    float interpolatedValue;

    std::vector<uint> vdb_prim_num;

    //vit_t.setCompressOnExit(true);

    for(vit.rewind(); !vit.atEnd(); vit.advanceTile() )
    {

      if( skipTile(vit, vdb_bbox, vdb_prim_num, voxel_size ,dst, 0 ) )  continue;

      //if(vdb_depth>0)  {  if( skipTile_vdb(vit, voxel_size, vdb_prim_num, vdb_grid, dst, vdb_depth+1)  ) continue;  }

      vit_t.setTile(vit);

      switch(operation)
      {
        case 0: // COPY
	  for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
	  {
	    interpolatedValue = 0;
	    dst->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );
	    const openvdb::Vec3R ijk(voxel_pos.x(), voxel_pos.y(), voxel_pos.z());

        for(j=0; j<vdb_prim_num.size(); j++)
	    {
        	interpolatedValue += sampler_vec[vdb_prim_num[j]].wsSample( ijk )*vdb_mult[ vdb_prim_num[j] ];
        }

        vit_t.setValue( interpolatedValue*src_mult );

            //vit_t.setValue( 1 );
	  }
	break;

        case 1:// ADD
	  for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
	  {
	    interpolatedValue = 0;
	    dst->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );
	    const openvdb::Vec3R ijk(voxel_pos.x(), voxel_pos.y(), voxel_pos.z());
	    for(j=0; j<vdb_prim_num.size(); j++)
	    {
	      interpolatedValue += sampler_vec[vdb_prim_num[j]].wsSample( ijk )*vdb_mult[ vdb_prim_num[j] ];
            }

	    vit_t.setValue( vit_t.getValue() + interpolatedValue*src_mult );
	  }
	break;

        case 2:// SUBTRACT
	  for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
	  {
	    interpolatedValue = 0;
	    dst->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );
	    const openvdb::Vec3R ijk(voxel_pos.x(), voxel_pos.y(), voxel_pos.z());
	    for(j=0; j<vdb_prim_num.size(); j++)
	    {
	      interpolatedValue += sampler_vec[vdb_prim_num[j]].wsSample( ijk )*vdb_mult[ vdb_prim_num[j] ];
	    }

	    vit_t.setValue( vit_t.getValue() - interpolatedValue*src_mult );
	  }
	break;

	case 3:
	  for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
	  {
	    interpolatedValue = 0;
	    dst->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );
	    const openvdb::Vec3R ijk(voxel_pos.x(), voxel_pos.y(), voxel_pos.z());

	    for(j=0; j<vdb_prim_num.size(); j++)
	    {
	      interpolatedValue += sampler_vec[vdb_prim_num[j]].wsSample( ijk )*vdb_mult[ vdb_prim_num[j] ];
	    }

	    vit_t.setValue( vit_t.getValue() * interpolatedValue*src_mult );
	  }
	break;

	case 4:
	  for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
	  {
	    interpolatedValue = 0;
	    dst->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );
	    const openvdb::Vec3R ijk(voxel_pos.x(), voxel_pos.y(), voxel_pos.z());

	    for(j=0; j<vdb_prim_num.size(); j++)
	    {
	      interpolatedValue += sampler_vec[vdb_prim_num[j]].wsSample( ijk )*vdb_mult[ vdb_prim_num[j] ];
	    }

	    vit_t.setValue( vit_t.getValue() / (interpolatedValue*src_mult) );
	  }
	break;

	case 5: // MAXIMUN
	  for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
	  {

//            if (vit_t.isStartOfTile() && boss->opInterrupt())
//                break;

	    interpolatedValue = -1000000;
	    dst->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );
	    const openvdb::Vec3R ijk(voxel_pos.x(), voxel_pos.y(), voxel_pos.z());


        //for(j=0; j<sampler_vec.size(); j++)
	    //{
	    //  interpolatedValue = std::max(  *vdbm, interpolatedValue ) ;
	    //}

		for(j=0; j<vdb_prim_num.size(); j++)
		{
			interpolatedValue = std::max( sampler_vec[ vdb_prim_num[j] ].wsSample( ijk )*src_mult*vdb_mult[ vdb_prim_num[j] ], interpolatedValue );
		}




//	    for(j=0; j<vdb_prim_num.size(); j++)
//	    {
//	      interpolatedValue = std::max(  sampler_vec[ vdb_prim_num[j] ].wsSample( ijk )*src_mult, interpolatedValue ) ;
//	    }
	    vit_t.setValue( std::max( vit_t.getValue() , interpolatedValue ) );
	  }
	break;

	// MINIMUM
	case 6:
	  for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
	  {
	    interpolatedValue = 100000;
	    dst->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );
	    const openvdb::Vec3R ijk(voxel_pos.x(), voxel_pos.y(), voxel_pos.z());

	    for(j=0; j<vdb_prim_num.size(); j++)
	    {
	      interpolatedValue = std::max( sampler_vec[ vdb_prim_num[j] ].wsSample( ijk )*src_mult*vdb_mult[ vdb_prim_num[j] ], interpolatedValue );
	    }
	    vit_t.setValue( std::min( vit_t.getValue() , interpolatedValue ) );
	  }
	break;

	// TEST

	case 7:
	  for (vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
	  {
	    vit_t.setValue( vdb_prim_num.size() );
	  }
	break;


      }

    }


  }


  // ---------------------------------------   CALCULATE VECTOR    ---------------------------------------------

template <class T1, class T2>  void m_gasvdbtofield<T1, T2>::calculateVectorPartial( SIM_VectorField* dst_x,
				                                 int idx,
				                                 std::vector< openvdb::VectorGrid::Ptr >& vdb_grid,
												 std::vector<UT_BoundingBoxF>& vdb_bbox,
				                                 std::vector<openvdb::math::Transform >& vdb_xform,
												 float src_mult,
												 int operation,
												 const UT_JobInfo &info )
  {

    UT_VoxelArrayIteratorF vit_x;
	UT_VoxelArrayIteratorF vit_y;
	UT_VoxelArrayIteratorF vit_z;

    UT_VoxelTileIteratorF  vit_t;

	std::vector<UT_Vector3F > vdb_vel(vdb_vel_in);

	SIM_RawField* velx = dst_x->getField(0);
	SIM_RawField* vely = dst_x->getField(1);
	SIM_RawField* velz = dst_x->getField(2);

	UT_VoxelArrayF* velx_rf = velx->fieldNC();
	UT_VoxelArrayF* vely_rf = vely->fieldNC();
	UT_VoxelArrayF* velz_rf = velz->fieldNC();

	vit_x.setArray(velx->fieldNC());
	vit_y.setArray(velx->fieldNC());
	vit_z.setArray(velx->fieldNC());

    UT_Vector3F voxel_pos;
    UT_Vector3F old_value;
    UT_Vector3F new_value;
    double old_val_length,new_val_length;


	UT_Vector3F source_pos;

    int i =0;


    int y_ix, y_iy, y_iz;
    int z_ix, z_iy, z_iz;

    UT_Vector3F voxel_size = velx->getVoxelSize();


    UT_Vector3F voxel_offset, voxel_offset_y, voxel_offset_z;


	if(dst_x->isFaceSampled())
		voxel_offset = voxel_size;
	else
		voxel_offset = UT_Vector3F(0,0,0);

	 voxel_offset_y = UT_Vector3F(voxel_offset.x()/2, -voxel_offset.y()/2 ,0);
	 voxel_offset_z = UT_Vector3F(voxel_offset.x()/2, 0 ,-voxel_offset.z()/2);

    std::vector<uint> vdb_prim_num;

    vit_x.setPartialRange(info.job(), info.numJobs());

    int j=0;
    openvdb::Vec3R val(0,0,0);
    openvdb::Vec3R tmp_vec(0,0,0);
    openvdb::Vec3R new_val(0,0,0);
    openvdb::Vec3R zero(0,0,0);

    UT_Vector3F velv(0,0,0);
	UT_Vector3F tmp_vec_h(0,0,0);
	UT_Vector3F h_zero(0,0,0);

	UT_Vector3F new_pos(0, 0, 0);
	double max_length;
	float  max_length2;

    //< openvdb::Vec3STree, openvdb::tools::PointSampler>
    std::vector< openvdb::tools::GridSampler< T1, T2> >  sampler_vec;

    for(i=0; i<vdb_xform.size(); i++)
    {
        openvdb::tools::GridSampler< T1, T2 >
        interpolator( vdb_grid[ i ]->constTree(), vdb_xform[i] );
        sampler_vec.push_back( interpolator );
    }


    for(vit_x.rewind(); !vit_x.atEnd(); vit_x.advanceTile() )
    {

		if( skipTile(vit_x, vdb_bbox, vdb_prim_num, voxel_size , velx, 0 ) )
			 continue;

		vit_t.setTile(vit_x);


		switch(operation)
		{
		case 0: //----------------------------------  VDB VALUES  ----------------------------------
			for(vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
		{
			velx->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );

			source_pos = voxel_pos + UT_Vector3F(0, -voxel_offset.y()/2 ,-voxel_offset.z()/2);

			dst_x->posToIndex( 1, voxel_pos + UT_Vector3F(voxel_offset.x()/2, -voxel_offset.y()/2 ,0),  y_ix, y_iy, y_iz );
			if( !vely_rf->isValidIndex(y_ix, y_iy, y_iz) ) continue;

			dst_x->posToIndex( 2, voxel_pos + UT_Vector3F(voxel_offset.x()/2, 0 ,-voxel_offset.z()/2), z_ix, z_iy, z_iz );
			if( !velz_rf->isValidIndex(z_ix, z_iy, z_iz) ) continue;

			old_value.assign(vit_t.getValue(), vely_rf->getValue(y_ix, y_iy, y_iz), velz_rf->getValue(z_ix, z_iy, z_iz) );

			const openvdb::Vec3R ijk(source_pos.x(), source_pos.y(), source_pos.z());

			val = zero;
			max_length  = 0.0;

			for(j=0; j<vdb_prim_num.size(); j++)
			{
				tmp_vec  = sampler_vec[vdb_prim_num[j]].wsSample( ijk )*vdb_mult[ vdb_prim_num[j] ];
				val += tmp_vec;
				max_length  = std::max( tmp_vec.length(),  max_length);
			}

			val.normalize();
			val *= max_length;

			new_value.assign(val[0], val[1], val[2]);

			new_value *= src_mult;

			new_val_length = new_value.length();
			old_val_length = old_value.length();

			if( (old_val_length + new_val_length) < vec_lenght_threshold )
			{
				//ADD
				new_value += old_value;
				vit_t.setValue( new_value.x() );
				vely_rf->setValue( y_ix, y_iy, y_iz, new_value.y() );
				velz_rf->setValue( z_ix, z_iy, z_iz, new_value.z() );

			}
			else
			{
				//MAXIMUM
				new_value += old_value;
				new_value.normalize();
				new_value *= std::max(new_val_length, old_val_length);
				vit_t.setValue( new_value.x() );
				vely_rf->setValue( y_ix, y_iy, y_iz, new_value.y() );
				velz_rf->setValue( z_ix, z_iy, z_iz, new_value.z() );
			}

		}
		break;

		case 1: //---------------------------------- V ATTRIBUTE  ----------------------------------

		for(vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
		{
			velx->indexToPos( vit_t.x(), vit_t.y(), vit_t.z(), voxel_pos  );

			source_pos = voxel_pos + UT_Vector3F(0, -voxel_offset.y()/2 ,-voxel_offset.z()/2);

			dst_x->posToIndex( 1, voxel_pos + UT_Vector3F(voxel_offset.x()/2, -voxel_offset.y()/2 ,0),  y_ix, y_iy, y_iz );
			dst_x->posToIndex( 2, voxel_pos + UT_Vector3F(voxel_offset.x()/2, 0 ,-voxel_offset.z()/2),  z_ix, z_iy, z_iz );

			if( !vely_rf->isValidIndex(y_ix, y_iy, y_iz) ) continue;
			if( !velz_rf->isValidIndex(z_ix, z_iy, z_iz) ) continue;

			old_value.assign(vit_t.getValue(), vely_rf->getValue(y_ix, y_iy, y_iz), velz_rf->getValue(z_ix, z_iy, z_iz) );

			const openvdb::Vec3R ijk(source_pos.x(), source_pos.y(), source_pos.z());

			//val         = zero;
			max_length  = 0.0;
			max_length2 = 0.0;
			velv        = h_zero;


			for(j=0; j<vdb_prim_num.size(); j++)
			{

				tmp_vec = sampler_vec[vdb_prim_num[j]].wsSample( ijk )*vdb_mult[ vdb_prim_num[j] ];

				tmp_vec_h = vdb_vel_in[vdb_prim_num[j]]*tmp_vec.x();

				velv += tmp_vec_h;

 			    max_length2 = std::max( tmp_vec_h.length(), max_length2);

			}


			velv.normalize();
			velv *= max_length2;

			new_value = velv;
			new_value *= src_mult;

			new_val_length = new_value.length();
			old_val_length = old_value.length();

			if( (old_val_length + new_val_length) < vec_lenght_threshold )
			{
				//ADD
				new_value += old_value;
				vit_t.setValue( new_value.x() );
				vely_rf->setValue( y_ix, y_iy, y_iz, new_value.y() );
				velz_rf->setValue( z_ix, z_iy, z_iz, new_value.z() );

			}
			else
			{
				//MAXIMUM
				new_value += old_value;
				new_value.normalize();
				new_value *= std::max(new_val_length, old_val_length);
				vit_t.setValue( new_value.x() );
				vely_rf->setValue( y_ix, y_iy, y_iz, new_value.y() );
				velz_rf->setValue( z_ix, z_iy, z_iz, new_value.z() );
			}

		}

		break;

		case 2: //----------------------------------  TEST  ----------------------------------


		for(vit_t.rewind(); !vit_t.atEnd(); vit_t.advance() )
		{
			vit_t.setValue(  vdb_prim_num.size() );
		}

		break;


		} //----------------------------------  END SWITCH  ---------------------------------


	}

  }


template <class T1, class T2>  void run_float_mt( SIM_RawField *dst,
				                                  std::vector< openvdb::FloatGrid::Ptr >& vdb_grid,
											      std::vector<UT_BoundingBoxF>& vdb_bbox,
				                                  std::vector<openvdb::math::Transform >& vdb_xform,
											      float src_mult,
											      int operation,
											      std::vector< float >& vdb_mult )
{
	m_gasvdbtofield<T1, T2>  ptf;
	ptf.vdb_mult  = vdb_mult;
	ptf.calculateScalar( dst, vdb_grid ,vdb_bbox, vdb_xform, src_mult, operation );
}



template <class T1, class T2>  void run_vector_mt( SIM_VectorField* dst_x,
				                                   int idx,
				                                   std::vector< openvdb::VectorGrid::Ptr >& vdb_grid,
				                                   std::vector< UT_Vector3F >  vdb_vel,
												   std::vector<UT_BoundingBoxF>& vdb_bbox,
				                                   std::vector<openvdb::math::Transform >& vdb_xform,
												   float src_mult,
												   int operation,
											       std::vector< float >& vdb_mult,
											       float vel_threshold )
{
	m_gasvdbtofield<T1, T2>  ptf;

	ptf.vdb_mult   = vdb_mult;
	ptf.vec_lenght_threshold = vel_threshold;
	ptf.vdb_vel_in = vdb_vel;;
    ptf.calculateVector( dst_x, 0, vdb_grid ,vdb_bbox, vdb_xform, src_mult, operation );

}



void
initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(vdb_to_field);
}


vdb_to_field::vdb_to_field(const SIM_DataFactory *factory)
    : BaseClass(factory)
{

}

vdb_to_field::~vdb_to_field()
{
}


const SIM_DopDescription *vdb_to_field::getDopDescription()
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

    static PRM_Name     timeName("time", "Time");
    static PRM_Default  timeDef(0,"$T");




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
        PRM_Template(PRM_ORD,    1, &vdbSamplerName,  PRMzeroDefaults, &vdbSamplerMenu),

        PRM_Template(PRM_FLT,    1, &srcMultName,   &srcMultDef, 0,   &srcMultRange),
        PRM_Template(PRM_FLT,    1, &velThresholdName,   &velThresholdDef, 0,   &velThresholdRange),
        PRM_Template(PRM_INT,    1, &substepsName,  &substepsDef, 0,  &substepsRange),
        PRM_Template(PRM_FLT_J,  1, &timeName,  &timeDef),
        PRM_Template()
    };

    static SIM_DopDescription    theDopDescription(
            true,                  // Should we make a DOP?
            "vdb_to_field",       // Internal name of the DOP.
            "vdb_to_field",       // Label of the DOP
            "$OS",              // Default data name
            classname(),           // The type of this DOP, usually the class.
            theTemplates);         // Template list for generating the DOP

    return &theDopDescription;
}

bool vdb_to_field::solveGasSubclass(SIM_Engine &engine,
									SIM_Object *obj,
									SIM_Time time,
									SIM_Time timestep)
{


	float src_mult = 0;

    if(getNormalize())
            src_mult = getSourceMultiply()* timestep * OPgetDirector()->getChannelManager()->getSamplesPerSec();
    else
            src_mult = getSourceMultiply();

    if(src_mult==0)	return true;

    SIM_ScalarField     *dstscalar;
    SIM_VectorField     *dstvector;
    SIM_ScalarField     *mask_field;

    SIM_DataArray   dst;
    SIM_DataArray   mask;

    getMatchingData(dst,  obj, GAS_NAME_FIELDDEST);
    getMatchingData(mask, obj, "source_mask");

    float stime = getTime();
    float gtime;
    int   susteps = getSubsteps();
    float _timestep = timestep;
    float vel_threshold =  getVelThreshold();
    int operation      =  getOperation();
    int voperation     =  getvOperation();



    for(int s=0; s<susteps; s++)
    {

	uint num_packed_vdb = 0;

    gtime = stime - (_timestep/susteps)*s;


    UT_String vdb_path;
    UT_String vdb_name;
    UT_String filed_name;
    GU_DetailHandle  packed_gdh;
    GU_Detail* packed_gdp;
    SOP_Node*        vdb_sop;
    OP_Context       vdb_context(gtime);

    getVdbPath(vdb_path);
    getVdbName(vdb_name);


    vdb_sop = getSOPNode( vdb_path ,1);
    if( !vdb_sop )
    {
      addError(obj, SIM_MESSAGE, "Invalid Vdb Path", UT_ERROR_WARNING);
      return true;
    }

    vdb_sop->forceRecook();

    packed_gdh = vdb_sop->getCookedGeoHandle(vdb_context);
    GU_DetailHandleAutoReadLock  objects_gdl(packed_gdh);
    packed_gdp = (GU_Detail*) objects_gdl.getGdp();

    if( !packed_gdp )
    {
      addError(obj, SIM_MESSAGE, "Invalid Vdb Source Path", UT_ERROR_WARNING);
      return true;
    }

    if(  !dst.entries() )
    {
      addError(obj, SIM_MESSAGE, "Fewer fields specified.", UT_ERROR_WARNING);
      return true;
    }



	dstscalar = SIM_DATA_CAST(dst(0), SIM_ScalarField);
	dstvector = SIM_DATA_CAST(dst(0), SIM_VectorField);

	std::vector< openvdb::FloatGrid::Ptr >  vdb_grig;
	std::vector< openvdb::VectorGrid::Ptr > vdb_grig_vector;
	std::vector<UT_BoundingBoxF>            vdb_bbox;
	std::vector<UT_BoundingBoxF>            vdb_nodes_bbox;
	std::vector< UT_Vector3F >              vdb_vel;
	std::vector< UT_Vector3F >              vdb_w;
	std::vector< UT_Vector3F >              vdb_center;
	std::vector< openvdb::Mat3R >           vdb_vel_xform;
	std::vector< float >                    vdb_mult;
	std::vector<openvdb::math::Transform >  vdb_xform;

    UT_Matrix4D mtrx;
    UT_Vector3F pt_pos;

	filed_name =  obj->getSubDataName(  obj->getSubDataIndex(dst(0)) );


    GEO_Primitive* prim;

    GU_ConstDetailHandle  vdb_gdh;
    GU_Detail*            vdb_gdp;
    GEO_Primitive*        vdb_gdp_prim;
	GEO_PrimVDB*          vbd_prim;

	GA_ROAttributeRef  nameAttrRef;
    GA_ROAttributeRef  multAttrRef;
	GA_ROAttributeRef  maxSubstepsAttrRef;

    std::vector< GEO_PrimVDB* > vdb_pntr_array;
    std::vector< UT_Matrix4D >  mtrx_array;
    std::vector< float >        mlt_array;
	std::vector< UT_Vector3F >  v_array;

    GA_ROAttributeRef mult_all_attr_ref, mult_name_attr_ref, v_attr_ref;

    GA_FOR_ALL_PRIMITIVES( packed_gdp, prim )
    {

		bool found_vdb = false;
		std::vector<int> prim_grids_num;

		maxSubstepsAttrRef = packed_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "max_substeps");
		if(maxSubstepsAttrRef)
		{
			GA_ROHandleI maxSubstepsAttr(maxSubstepsAttrRef.getAttribute());
			int max_substeps = maxSubstepsAttr.get(prim->getMapOffset());
			int sd = susteps/max_substeps;

			//std::cout << "max_substeps: " << max_substeps << " susteps: " << susteps << " sd: " <<  sd <<  " s: " << s << "  " << "\n";

			if( (max_substeps<susteps) && (max_substeps!=0) && (sd!=0) )
				if( s%sd != 0 ) continue;
		}

    	if(prim->getTypeId() == GA_PRIMINTERNALSENTINEL)
    	{

    		GU_PrimPacked* packed_prim = UTverify_cast<GU_PrimPacked*>( prim );

			vdb_gdh = packed_prim->getPackedDetail();

			GU_DetailHandleAutoReadLock  objects_gdl(vdb_gdh);
            vdb_gdp = (GU_Detail*) objects_gdl.getGdp();

			nameAttrRef = vdb_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "name");
		    if( nameAttrRef.isInvalid() )
		    {
				UT_String msg;
                msg.sprintf("\nPacked primitive #%0d doesn't contain inside primitive attribute \"name\" \n", prim->getMapOffset() );

				addError(obj, SIM_MESSAGE, msg, UT_ERROR_WARNING);
				continue;
		    }
			GA_ROHandleS nameAttr(nameAttrRef.getAttribute());
			

			GA_FOR_ALL_PRIMITIVES( vdb_gdp, vdb_gdp_prim )
			{

				if( vdb_gdp_prim->getTypeId() != GA_PRIMVDB ) continue;
				if( vdb_name != nameAttr.get(vdb_gdp_prim->getMapOffset()) ) continue;

				GEO_PrimVDB* vbd_prim = UTverify_cast<GEO_PrimVDB*>( vdb_gdp_prim );

				if(  dstscalar && vbd_prim->getStorageType() != UT_VDB_FLOAT ) continue;
				if(  dstvector && vbd_prim->getStorageType() != UT_VDB_VEC3F ) continue;

				num_packed_vdb++;

				vdb_pntr_array.push_back(vbd_prim);

				packed_prim->getFullTransform4(mtrx);
				mtrx_array.push_back(mtrx);

				mult_all_attr_ref  = packed_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "multiply");
				v_attr_ref         = packed_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "v");
				mult_name_attr_ref = packed_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, vdb_name);

				float total_mult = 1.0;

				if( mult_all_attr_ref.isValid() )
				{
					GA_ROHandleF mult_attr(mult_all_attr_ref);
                    total_mult *= mult_attr.get(prim->getMapOffset());
				}

				if( mult_name_attr_ref.isValid() )
				{
					GA_ROHandleF mult_attr(mult_name_attr_ref);
                    total_mult *= mult_attr.get(prim->getMapOffset());
				}

				if( v_attr_ref.isValid() )
				{
					GA_ROHandleV3 v_attr(v_attr_ref);
                    v_array.push_back( v_attr.get(prim->getMapOffset()) );
				}
				else
				{
					v_array.push_back( UT_Vector3F(0,0,0) );
				}

				mlt_array.push_back(total_mult);
			}

    	}
    	else if(  prim->getTypeId() == GA_PRIMVDB )
    	{
    		nameAttrRef = packed_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "name");
    		if( nameAttrRef.isInvalid() )
		    {
				UT_String msg;
                msg.sprintf("name Attribute on VDB doesn't exist", prim->getMapOffset() );

				addError(obj, SIM_MESSAGE, msg, UT_ERROR_WARNING);
				continue;
		    }


			GA_ROHandleS nameAttr(nameAttrRef.getAttribute());
			if( vdb_name != nameAttr.get(prim->getMapOffset()) ) continue;

			GEO_PrimVDB* vbd_prim = UTverify_cast<GEO_PrimVDB*>( prim );

			if(  dstscalar && vbd_prim->getStorageType() != UT_VDB_FLOAT ) continue;
			if(  dstvector && vbd_prim->getStorageType() != UT_VDB_VEC3F ) continue;

			num_packed_vdb++;

			vdb_pntr_array.push_back(vbd_prim);

			UT_Matrix4D mtrx;
			mtrx.identity();

			mtrx_array.push_back(mtrx);

			mult_all_attr_ref  = packed_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "multiply");
			v_attr_ref         = packed_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "v");
			mult_name_attr_ref = packed_gdp->findAttribute(GA_ATTRIB_PRIMITIVE, vdb_name);

			float total_mult = 1.0;

			if( mult_all_attr_ref.isValid() )
			{
				GA_ROHandleF mult_attr(mult_all_attr_ref);
                total_mult *= mult_attr.get(prim->getMapOffset());
			}

			if( mult_name_attr_ref.isValid() )
			{
				GA_ROHandleF mult_attr(mult_name_attr_ref);
                total_mult *= mult_attr.get(prim->getMapOffset());
			}

			if( v_attr_ref.isValid() )
			{
				GA_ROHandleV3 v_attr(v_attr_ref);
                v_array.push_back( v_attr.get(prim->getMapOffset()) );
			}
			else
			{
				v_array.push_back( UT_Vector3F(0,0,0) );
			}

			mlt_array.push_back(total_mult);

    	}


	}


	if(dstscalar)
	{

		for(int i=0; i<vdb_pntr_array.size(); i++)
		{

		if( mlt_array[i] == 0.0 ) continue;

		GEO_PrimVDB* vbd_prim = vdb_pntr_array[i];

        openvdb::GridBase::Ptr vdb_base_grid = vbd_prim->getGridPtr();
        openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(vdb_base_grid);

        if( grid.get()->empty() ) continue;

        vdb_grig.push_back( grid );


        openvdb::math::CoordBBox bb;
        grid.get()->tree().evalLeafBoundingBox(bb);

        openvdb::CoordBBox box = grid->evalActiveVoxelBoundingBox();
        openvdb::Coord minmaxIndex[2]={box.min(),box.max()};
        openvdb::Vec3d minmaxWorld[2]={openvdb::Vec3f(FLT_MAX,FLT_MAX,FLT_MAX),openvdb::Vec3f(-FLT_MAX,-FLT_MAX,-FLT_MAX)};


        for(int i=0;i<8;i++)
        {
            openvdb::Vec3f v(minmaxIndex[i&1].x(),minmaxIndex[(i&2)>>1].y(),minmaxIndex[(i&4)>>2].z());
            openvdb::math::Vec3d vworld=grid->indexToWorld(v);

            for(int k=0;k<3;k++)
            {
                    minmaxWorld[0][k]=std::min(minmaxWorld[0][k],vworld[k]);
                    minmaxWorld[1][k]=std::max(minmaxWorld[1][k],vworld[k]);
            }
        }

        float mult;

		mtrx = mtrx_array[i];
		openvdb::math::Mat4d mat(mtrx.data());

		if( !openvdb::math::isAffine(mat) || !openvdb::math::isInvertible(mat) )
		{
		    UT_String msg;
            msg.sprintf("Matrix on point #%0d in not Affine or not Invertible", -1 );

            addError(obj, SIM_MESSAGE, msg, UT_ERROR_WARNING);

            mtrx.identity();
            mtrx.setTranslates(UT_Vector3F(i,0,0));
            mat = openvdb::math::Mat4d(mtrx.data());
		}

		openvdb::math::Transform xform = grid->transform();
		xform.postMult(mat);
		vdb_xform.push_back( xform );

        UT_BoundingBoxF bbox(   UT_Vector3F( minmaxWorld[0].x(), minmaxWorld[0].y(), minmaxWorld[0].z() ) ,
                                UT_Vector3F( minmaxWorld[1].x(), minmaxWorld[1].y(), minmaxWorld[1].z() ) );


        bbox.transform(mtrx);

        vdb_bbox.push_back( bbox );
        vdb_mult.push_back( mlt_array[i] );

        }

	}


	if(dstvector)
	{
		for(int i=0; i<vdb_pntr_array.size(); i++)
		{

			if( mlt_array[i] == 0.0 ) continue;

			GEO_PrimVDB* vbd_prim = vdb_pntr_array[i];

			openvdb::GridBase::Ptr vdb_base_grid =  vbd_prim->getGridPtr();
            openvdb::VectorGrid::Ptr grid = openvdb::gridPtrCast<openvdb::VectorGrid>(vdb_base_grid);

		    if( grid.get()->empty() )  continue;

			openvdb::Mat3R mtrx_ident;
			mtrx_ident.setIdentity();

			vdb_grig_vector.push_back( grid );

            openvdb::math::CoordBBox bb;
            grid.get()->tree().evalLeafBoundingBox(bb);

            openvdb::CoordBBox box = grid->evalActiveVoxelBoundingBox();
            openvdb::Coord minmaxIndex[2]={box.min(),box.max()};
            openvdb::Vec3d minmaxWorld[2]={openvdb::Vec3f(FLT_MAX,FLT_MAX,FLT_MAX),openvdb::Vec3f(-FLT_MAX,-FLT_MAX,-FLT_MAX)};
            for(int i=0;i<8;i++)
            {
                openvdb::Vec3f v(minmaxIndex[i&1].x(),minmaxIndex[(i&2)>>1].y(),minmaxIndex[(i&4)>>2].z());
                openvdb::math::Vec3d vworld=grid->indexToWorld(v);
                for(int k=0;k<3;k++)
                {
                    minmaxWorld[0][k]=std::min(minmaxWorld[0][k],vworld[k]);
                    minmaxWorld[1][k]=std::max(minmaxWorld[1][k],vworld[k]);
                }
            }


			vdb_mult.push_back( mlt_array[i] );
			vdb_vel.push_back(  v_array[i] );


            mtrx = mtrx_array[i];
            openvdb::math::Mat4d mat(mtrx.data());

            if(  !openvdb::math::isAffine(mat) || !openvdb::math::isInvertible(mat) )
            {
                UT_String msg;
                msg.sprintf("Matrix on point #%0d in not Affine or not Invertible", -1 );

                addError(obj, SIM_MESSAGE, msg, UT_ERROR_WARNING);

                mtrx.identity();
                mtrx.setTranslates(pt_pos);
                mat = openvdb::math::Mat4d(mtrx.data());
            }

            openvdb::math::Transform xform = grid->transform();
            xform.postMult(mat);

            vdb_xform.push_back( xform );

            UT_BoundingBoxF bbox(   UT_Vector3F( minmaxWorld[0].x(), minmaxWorld[0].y(), minmaxWorld[0].z() ) ,
                                    UT_Vector3F( minmaxWorld[1].x(), minmaxWorld[1].y(), minmaxWorld[1].z() ) );

            bbox.transform(mtrx);
            vdb_bbox.push_back( bbox );

		}

	}

	 UT_String msg;
     msg.sprintf("\n\nUsing %0d vdb out of %0d", vdb_xform.size(), num_packed_vdb);

	addError( obj, SIM_MESSAGE, msg, UT_ERROR_MESSAGE );

    if(dstscalar && vdb_grig.size()>0 )
    {
      	  SIM_RawField *dstfield = dstscalar->getField();

      	  switch( getVdbSampler() )
      	  {
      	  	case 0:
      	  		run_float_mt<openvdb::FloatTree, openvdb::tools::PointSampler>( dstfield, vdb_grig ,vdb_bbox, vdb_xform, src_mult, operation, vdb_mult );
      	  	break;
      	  	case 1:
      	  		run_float_mt<openvdb::FloatTree, openvdb::tools::BoxSampler>( dstfield, vdb_grig, vdb_bbox, vdb_xform, src_mult, operation, vdb_mult );
      	  	break;
      	  	case 2:
      	   		run_float_mt<openvdb::FloatTree, openvdb::tools::QuadraticSampler>( dstfield, vdb_grig ,vdb_bbox, vdb_xform, src_mult, operation, vdb_mult );
      	    break;

      	  }

      	  vdb_grig.clear();
		  vdb_bbox.clear();
		  vdb_xform.clear();
		  vdb_mult.clear();
    }

    if(dstvector && vdb_grig_vector.size()>0)
    {
    	switch( getVdbSampler() )
      	{
      		case 0:
				run_vector_mt<openvdb::VectorGrid, openvdb::tools::StaggeredPointSampler>( dstvector, 0,  vdb_grig_vector, vdb_vel, vdb_bbox, vdb_xform, src_mult, voperation, vdb_mult, vel_threshold  );
      		break;
      		case 1:
				run_vector_mt<openvdb::VectorGrid, openvdb::tools::StaggeredBoxSampler>( dstvector, 0, vdb_grig_vector, vdb_vel, vdb_bbox, vdb_xform, src_mult, voperation, vdb_mult, vel_threshold  );
      		break;
      		case 2:
				run_vector_mt<openvdb::VectorGrid, openvdb::tools::StaggeredQuadraticSampler>( dstvector, 0,  vdb_grig_vector, vdb_vel,vdb_bbox, vdb_xform, src_mult, voperation, vdb_mult, vel_threshold  );
      		break;

      	}

      	vdb_grig_vector.clear();
      	vdb_bbox.clear();
		vdb_xform.clear();
		vdb_mult.clear();

    }

       if (dstscalar)
	dstscalar->pubHandleModification();
       if (dstvector)
	dstvector->pubHandleModification();

    }

    return true;
}


