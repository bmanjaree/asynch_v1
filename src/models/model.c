#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>

#include <models/model.h>


/// Models
extern const AsynchModel model_252;
extern const AsynchModel model_190;


AsynchModel const * GetModel(unsigned short model_uid)
{
    AsynchModel const *res = NULL;

    //// An array of available model definitions
    //AsynchModel const * models[] = {
    //    &model_190,
    //    &model_252
    //};

    ///// Number of model registered
    //unsigned int num_models = 0;

    //// Search for model
    //for (unsigned int i = 0; i < num_models; i++)
    //{
    //    if (models[i]->uid == model_uid)
    //    {
    //        res = models[i];
    //        break;
    //    }
    //}

    return res;
}
