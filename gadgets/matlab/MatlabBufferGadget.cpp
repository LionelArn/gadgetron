#include "MatlabBufferGadget.h"
#include "MatlabUtils.h"

std::mutex mutex_;

namespace Gadgetron{

int MatlabBufferGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
{
    GDEBUG("Starting MatlabBufferGadget::process");
    
    std::lock_guard<std::mutex> lock(mutex_);   

	// Initialize a string for matlab commands
	std::string cmd;

	auto recon_data = m1->getObjectPtr();
	mwSize nencoding_spaces = recon_data->rbit_.size();

	const char* fieldnames[2] = {"data","reference"};
	auto reconArray = mxCreateStructArray(1,&nencoding_spaces,2,fieldnames); // what is this mysterious encoding_spaces ?

    // 2e9 bytes data is the published (as of 2017a) hardcoded limit that engPutVariable can transfer.
    // Empirically, it seems that variables up to 2^32 bytes (~4.3 GB) can be sent.
    size_t max_data_size = 2e9;
    
    // compute the data size in bytes (omitting the reference part)
    // There is probably a cleaner way to do it.
    // get_number_of_bytes()
    size_t data_bytes = 0;
    for (int i=0; i < recon_data->rbit_.size(); i++)
    {
        size_t bytes = sizeof(recon_data->rbit_[i].data_.data_[0]);   
        for(int j=0; j<recon_data->rbit_[i].data_.data_.get_number_of_dimensions(); ++j)
            bytes *= recon_data->rbit_[i].data_.data_.get_size(j);
        data_bytes += bytes;
    }
    
    GDEBUG("Bucket size: %lu bytes\n", data_bytes);
    
    if(0)//data_bytes < max_data_size) 
    {
        // the dataset is small enough to be sent all at once (original code)
        for (int i = 0; i <  recon_data->rbit_.size(); i++)
        {
            auto mxrecon = BufferToMatlabStruct(&recon_data->rbit_[i].data_);
            mxSetField(reconArray,i,"data",mxrecon);
            if (recon_data->rbit_[i].ref_)
            {
                auto mxref = BufferToMatlabStruct(recon_data->rbit_[i].ref_.get_ptr());
                mxSetField(reconArray,i,"reference",mxref);
            }
        }
        
        GDEBUG("Sending the whole buckets...\n");
        engPutVariable(engine_, "recon_data", reconArray);
        GDEBUG("done\n");
    }
    else
    {
        // the dataset needs to be sent in multiple packets
        // The algorithm here splits the multidimensional arrays (data.data
        // and reference.data) into n_packets in the RO dimension. After all
        // packets are sent, MATLAB reconcatenates everything.
        
        int n_packets = 2;//ceil( float(data_bytes) / float(max_data_size) );
        
        GDEBUG("Bucket size limit reached, parsing it into %i packets.\n", n_packets);
        
        for (int i = 0; i <  recon_data->rbit_.size(); i++)
        {            
            // Create the regular MATLAB structure, but omits the data for the fields "data" and "reference".
            auto mxrecon = BufferToMatlabStruct(&recon_data->rbit_[i].data_, true);
            mxSetField(reconArray,i,"data",mxrecon);
            if (recon_data->rbit_[i].ref_)
            {
                auto mxref = BufferToMatlabStruct(recon_data->rbit_[i].ref_.get_ptr(), true);
                mxSetField(reconArray,i,"reference",mxref);
            }
            
            
            // send the packets
            //size_t n_RO = sizeof(recon_data->rbit_[i].data_.data_) / sizeof(recon_data->rbit_[i].data_.data_[0]);
            
            float step = float(recon_data->rbit_[i].data_.data_.get_size(0))/float(n_packets);
            
            GDEBUG("Starting to process packets for index %i:\n", i+1);
            for(int p = 0; p < n_packets; p++)
            {

                
                // create the packet. A copy of the data is being done here,
                // which overall increase the RAM usage if packets are needed.
                // There may be a more efficient way to do this.
                
                // (RO) indexes of data to be split
                size_t beg = roundf(float(p  )*step       );
                size_t end = roundf(float(p+1)*step - 1.0f);
                GDEBUG("Creating data packet #%i: from index %lu to %lu...\n", p+1, (long unsigned) beg, (long unsigned) end);
                
                size_t dim_1_n_elem =   recon_data->rbit_[i].data_.data_.get_size(1)*
                                        recon_data->rbit_[i].data_.data_.get_size(2)*
                                        recon_data->rbit_[i].data_.data_.get_size(3)*
                                        recon_data->rbit_[i].data_.data_.get_size(4)*
                                        recon_data->rbit_[i].data_.data_.get_size(5)*
                                        recon_data->rbit_[i].data_.data_.get_size(6);
                
                //GDEBUG("Dim 1 element bytes: %lu\n", (long unsigned) bytes_dim_1);
                /*
                std::complex<float> packet  [end-beg]
                                            [recon_data->rbit_[i].data_.data_.get_size(1)]
                                            [recon_data->rbit_[i].data_.data_.get_size(2)]
                                            [recon_data->rbit_[i].data_.data_.get_size(3)]
                                            [recon_data->rbit_[i].data_.data_.get_size(4)]
                                            [recon_data->rbit_[i].data_.data_.get_size(5)]
                                            [recon_data->rbit_[i].data_.data_.get_size(6)];
                 std::copy(  &(recon_data->rbit_[i].data_.data_[0]) + beg*bytes_dim_1,
                             &(recon_data->rbit_[i].data_.data_[0]) + end*bytes_dim_1,
                             (packet[0][0][0][0][0][0][0]));
                 */
                
                /*
                mwSize ndim = input->get_number_of_dimensions();
                mwSize* dims = new mwSize[ndim];
                for (size_t i = 0; i < ndim; i++)
                    dims[i] = input->get_size(i);
                T* raw_data = (T*) mxCalloc(input->get_number_of_elements(),sizeof(T));
                memcpy(raw_data,input->get_data_ptr(),input->get_number_of_bytes());
                auto result =  mxCreateNumericMatrix(0,0,MatlabClassID<T>::value,isComplex<T>::value);
                mxSetDimensions(result,dims,ndim);
                mxSetData(result,raw_data);
                return result;
                */
                
                // for some reason, this methods only sends some weird real data to matlab.
                /*
                std::complex<float> * packet = (std::complex<float> *) mxCalloc((end-beg+1)*(bytes_dim_1/sizeof(std::complex<float>)), sizeof(std::complex<float>));
                memcpy(&(packet[0]), &(recon_data->rbit_[i].data_.data_[beg]),  (end-beg+1)*bytes_dim_1);
                
                // convert the packet to MALTAB array
                mwSize packet_ndim = recon_data->rbit_[i].data_.data_.get_number_of_dimensions();
                mwSize* packet_dims = new mwSize[packet_ndim];
                packet_dims[0] = end-beg+1;
                for (size_t j = 1; j < packet_ndim; j++)
                    packet_dims[j] = recon_data->rbit_[i].data_.data_.get_size(j);

                auto mxdata =  mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxCOMPLEX);
                mxSetDimensions(mxdata,packet_dims,packet_ndim);
                mxSetData(mxdata,packet);
                */
                
                
                /*
                size_t ndim = input->get_number_of_dimensions();
                mwSize* dims = new mwSize[ndim];
                for (size_t i = 0; i < ndim; i++)
                    dims[i] = input->get_size(i);

                REAL* real_data = (REAL*) mxCalloc(input->get_number_of_elements(),sizeof(REAL));
                REAL* imag_data = (REAL*) mxCalloc(input->get_number_of_elements(),sizeof(REAL));

                complext<REAL>* raw_data = input->get_data_ptr();
                for (size_t i = 0; i < input->get_number_of_elements(); i++){
                    real_data[i] = real(raw_data[i]);
                    imag_data[i] = imag(raw_data[i]);
                }

                auto result  =  mxCreateNumericMatrix(0,0,MatlabClassID<REAL>::value,isComplex<complext<REAL>>::value);
                mxSetDimensions(result,dims,ndim);
                mxSetData(result,real_data);
                mxSetImagData(result,imag_data);

                auto ndims_test = mxGetNumberOfDimensions(result);
                */
                
                size_t packet_n_elem = (end-beg+1) * dim_1_n_elem;
                
                size_t packet_ndim = recon_data->rbit_[i].data_.data_.get_number_of_dimensions();
                mwSize* packet_dims = new mwSize[packet_ndim];
                packet_dims[0] = end-beg+1;
                for (size_t j = 1; j < packet_ndim; j++)
                    packet_dims[j] = recon_data->rbit_[i].data_.data_.get_size(j);

                float* real_data = (float*) mxCalloc(packet_n_elem, sizeof(float));
                float* imag_data = (float*) mxCalloc(packet_n_elem, sizeof(float));
                
                size_t start = beg*dim_1_n_elem;
                
                GDEBUG("Index: start %lu, end: %lu\n", start, start + packet_n_elem - 1);
                
                std::complex<float>* raw_data = recon_data->rbit_[i].data_.data_.get_data_ptr();
                for (size_t j = 0; j < packet_n_elem; j++){
                    real_data[j] = real(raw_data[start + j]);
                    imag_data[j] = imag(raw_data[start + j]);
                }
                
                //for (size_t j = 0; j < start + packet_n_elem -1 + 2e9; j+=1000){
                //    GDEBUG("index: %lu: %f + %f*i\n", j, real(recon_data->rbit_[i].data_.data_[j]), imag(recon_data->rbit_[i].data_.data_[j]));
                //}
                    
                auto mxdata =  mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxCOMPLEX);
                mxSetDimensions(mxdata, packet_dims, packet_ndim);
                mxSetData      (mxdata, real_data);
                mxSetImagData  (mxdata, imag_data);
                
                GDEBUG("Sending data packet #%i...\n", p+1);
                
                std::string cmd = "data_" + std::to_string(i) + "_" + std::to_string(p);
                engPutVariable(engine_, cmd.c_str(), mxdata);
                
                std::string wcmd = "fprintf(2,evalc('whos'))";
                send_matlab_command(wcmd);
                
                wcmd = "figure; imagesc(abs(squeeze(" + cmd + "(size(" + cmd + ",1)/2,:,:,1)))); drawnow;"+
                       "figure; imagesc(abs(squeeze(" + cmd + "(:,128,:,1)))); drawnow;"+
                       "figure; imagesc(abs(squeeze(" + cmd + "(:,:,128,1)))); drawnow; pause(10);";
                send_matlab_command(wcmd);
                
                /*
                // do the same for the reference
                if (recon_data->rbit_[i].ref_)
                {
                    GDEBUG("Creating reference packet #%i...\n", p+1);
                    // create the ref packet
                    auto packet_ref = malloc( (end-beg)*sizeof(recon_data->rbit_[i].ref_.data_[0]) );
                    std::copy( &(recon_data->rbit_[i].ref_.data_[beg]),
                               &(recon_data->rbit_[i].ref_.data_[end]),
                               &(packet_ref[0]));
                    auto mxdata_ref = hoNDArrayToMatlab(&packet_ref);

                    // convert the ref packet to MALTAB array and send it
                    GDEBUG("Sending reference packet #%i...\n", p+1);
                    engPutVariable(engine_, "ref_" + std::to_string(i) + "_" + std::to_string(p), mxdata_ref);
                    free(packet_ref);
                }*/
            }
        }
        engPutVariable(engine_, "recon_data", reconArray);
        
        /*
        //send the command to reconcatenate the data and ref
        for (int i = 0; i <  recon_data->rbit_.size(); i++)
        {
            GDEBUG("MATLAB concatenation for index %i...\n", i+1);
            
            // create a concatenation MATLAB command
            string concat_data = "[";
            string concat_ref  = "[";
            for(int p = 0; p < n_packets; p++)
            {
                concat_data += "data_" + std::to_string(i) + "_" + std::to_string(p) + "; ";
                if (recon_data->rbit_[i].ref_)
                    concat_data += "ref_" + std::to_string(i) + "_" + std::to_string(p) + "; ";
            }
            
            // send the concatenation command to MATLAB
            send_matlab_command("recon_data.data(" + std::to_string(i) + ").data  = " + concat_data + "];");
            if (recon_data->rbit_[i].ref_)
                send_matlab_command("recon_data.ref(" + std::to_string(i) + ").data  = " + concat_ref + "];");
            
            // clear the MATLAB data copies
            for(int p = 0; p < n_packets; p++)
            {
                send_matlab_command("clear " + "data_" + std::to_string(i) + "_" + std::to_string(p) + "; ");
                if (recon_data->rbit_[i].ref_)
                    send_matlab_command("clear " + "ref_" + std::to_string(i) + "_" + std::to_string(p) + "; ");
            }
        }        */
    }
    
    cmd = "[imageQ,bufferQ] = matgadget.run_process(recon_data); matgadget.emptyQ();";
    send_matlab_command(cmd);


	// Get the size of the gadget's queue
	mxArray *imageQ = engGetVariable(engine_, "imageQ");
	if (imageQ == NULL) {
		GERROR("Failed to get the imageQ from matgadget\n");
		return GADGET_FAIL;
	}

	size_t qlen = mxGetNumberOfElements(imageQ);
	if (debug_mode_) {
	GDEBUG("Image Queue size: %d \n", qlen);
	}
	const mwSize* dims = mxGetDimensions(imageQ);
	mwSize ndims = mxGetNumberOfDimensions(imageQ);

	GDEBUG("Number of ndims %i \n",ndims);

	//Read all Image bytes
	for (mwIndex idx = 0; idx < qlen; idx++) {
		mxArray *res_hdr  = mxGetField(imageQ, idx, "bytes");
		mxArray *res_data = mxGetField(imageQ, idx, "image");

		GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 =
				new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		ISMRMRD::ImageHeader *hdr_new = m3->getObjectPtr();
		memcpy(hdr_new, mxGetData(res_hdr), sizeof(ISMRMRD::ImageHeader));

		auto image= MatlabToHoNDArray<std::complex<float>>(res_data);
		auto m4 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >(image);
		auto dims = *image.get_dimensions();

		m3->cont(m4);
		if (this->next()->putq(m3) < 0) {
			GDEBUG("Failed to put Image message on queue\n");
			return GADGET_FAIL;
		}

	}
	//Match engGetVariable with mxDestroy___s


	mxArray* bufferQ = engGetVariable(engine_,"bufferQ");

	qlen = mxGetNumberOfElements(bufferQ);
	if (debug_mode_) {
		GDEBUG("Buffer Queue size: %d \n", qlen);
		}

	for (mwIndex idx = 0; idx <qlen; idx++){

		IsmrmrdReconData output_data;
		IsmrmrdReconBit bit;
		bit.data_ = MatlabStructToBuffer(mxGetField(bufferQ,idx,"data"));

		auto ref = mxGetField(bufferQ,idx,"reference");
		if (ref){
			GDEBUG("Adding reference");
			bit.ref_ = MatlabStructToBuffer(ref);
		}
		output_data.rbit_.push_back(bit);
		auto m3 = new GadgetContainerMessage<IsmrmrdReconData>(output_data);
		if (this->next()->putq(m3) < 0){
			GDEBUG("Failed to put Buffer message on queue\n");
			return GADGET_FAIL;
		}

	}





	mxDestroyArray(bufferQ);
	mxDestroyArray(imageQ);
	mxDestroyArray(reconArray); //We're not supposed to delete this?

	// We are finished with the incoming messages m1 and m2
	m1->release();

	return GADGET_OK;
}
}


