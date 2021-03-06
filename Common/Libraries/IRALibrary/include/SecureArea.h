#ifndef _SECUREAREA_H
#define _SECUREAREA_H

/* **************************************************************************************************** */
/* IRA Istituto di Radioastronomia                                                                      */
/*                                                                                                      */
/* This code is under GNU General Public Licence (GPL).                                                 */
/*                                                                                                      */
/* Who                                when            What                                              */
/* Andrea Orlati(aorlati@ira.cnr.it)  03/01/2005      Creation                                          */
/* Andrea Orlati(aorlati@ira.cnr.it)  22/02/2005      Restructured in order to avoid deadlock if a mutex*/ 
/*                                                    was not released     								*/
/* Andrea Orlati(aorlati@ira.cnr.it)  16/08/2005      Added Init method to CSecureArea					*/ 


#include <new>
#include "baciThread.h"

namespace IRA {

template <class >class CSecureArea;

/**
 * This template class is used by the <i>CSecureArea</i> to make available the resource protected by a mutex mechanism.
 * This class is not intended to be used "stand-alone" and has been introduced in order to avoid deadlock if the mutex was not 
 * released. An object CSecureArea can be used as it is a pointer to the protected object (template typename). If the mutex is not 
 * released esplicitally it will be released when the <i>CSecAreaResourceWrapper</i> object get out of scope. For more information
 * see the <i>CSecureArea</i> documentation.
 * @author <a href=mailto:a.orlati@ira.cnr.it>Andrea Orlati</a>,
 * Istituto di Radioastronomia, Italia
 * <br>
*/	
template <class X> class CSecAreaResourceWrapper {
public:
	friend class CSecureArea<X>;
	/**
	 * Copy constructor
	 */
	CSecAreaResourceWrapper(const CSecAreaResourceWrapper& rSrc);
   /**
	 * Destructor
   */
	~CSecAreaResourceWrapper();
	/**
	 * Releases esplicitally the lock if it is still acquired, otherwise it does nothing.
	*/
	void Release();

	/**
	 * @return the pointer of the protected resource. If the lock is already released that means that
	 *         th program has a bug and the process is terminated (assert)
	*/
	X* operator->();
	/**
	 * @return the reference to the protected resource. If the lock is already released that means that
	 *         the program has a bug and the process is terminated (assert)
	*/
	X& operator*();

	/**
	 * it could be used to get direct access to the protected resource....very dangerous because
	 * you could access a resource which is already freed
	 */
	operator const X *() const {
		assert(m_pLockResource!=NULL);
		return m_pLockResource;
	}
private:
	/** Pointer to the mutex */
	BACIMutex* m_pLockMutex;
	/** Pointer to the protected resource */
	X* m_pLockResource;
   /** Indicates if the lock has been released or not */
	bool m_bReleased;
	/**
	 * Constructor.
	 * @param Mutex pointer to the mutex that must be allocated in order to synchronize the resource.
	 * @param Res pointer to the resource
	*/
	CSecAreaResourceWrapper(BACIMutex *Mutex,X *Res);
	/**
	 * Copy operator
	 */
	CSecAreaResourceWrapper<X>& operator=(const CSecAreaResourceWrapper& rSrc); // no implementation
};

/**
 * This template class contains the mechanism that can be used to protect a resource and synchronize all the threads that need to
 * make use of it. A thread must call <i>Get()</i> to be granted the priviledge to use the resource. An object <i>CSecAreaResourceWrapper</i>
 * is returned that let the developer use the resources and eventually forget to release the mutex lock. An example of usage:
 * <pre>
 *		//allocate the secure area and protect an integer.....
 *		CSecureArea<int>Area(new int(10));
 *		//obtain the access to the resource....
 *		CSecAreaResourceWrapper<int> Value=Area.Get();
 *		//make use of it....
 *		printf("%d",*Value);
 *		*Value=40;
 * </pre>
 * @author <a href=mailto:a.orlati@ira.cnr.it>Andrea Orlati</a>,
 * Istituto di Radioastronomia, Italia
 * <br>
*/
template <class T> class CSecureArea {
public:
	/**
	 * Constructor
	 * @param Obj pointer to the object/resource that must be protected
   */
	CSecureArea(T* Obj);
	/** Constructor
	 * @param Alloc this boolean tells the constructor if it has to allocate the resource automatically or not.
	*/
	CSecureArea(bool Alloc=false);
	/**
    * Destructor.
    * It also takes charge of deleting the resource.
	*/
	~CSecureArea();
	/**
    * Call this member function to pass the pointer to the object you need to protect. If the object/resoruce was already 
	* initialized, this function has no effect.
	* @param Obj pointer to the object/resource that must be protected
	*/
	void Init(T* Obj);
	/**
    * Call this member function to be granted the possibility to access the resource. If the lock of the resource is already
    * allocated this function will wait until it is released.
    * @return the resource wrapper that permits the developer to make use of the property.
	*/
	CSecAreaResourceWrapper<T> Get();

	/**
	* This method has been added for debugging purposes. The semantics is the same of the standard <i>Get</i> method but it also prints some information about the mutex and its owner.
    * Call this member function to be granted the possibility to access the resource. If the lock of the resource is already
    * allocated this function will wait until it is released.
    * @param owner string identifier of the mutex owner
    * @return the resource wrapper that permits the developer to make use of the property.
	*/
	CSecAreaResourceWrapper<T>Get(IRA::CString owner,bool show=false);

private:
	/** Mutex used to lock/unlock the resource */
	BACIMutex m_Mutex;
   	/** Pointer to the resource */
	T *m_pResource;
	IRA::CString m_resourceOwner;
	CSecureArea(const CSecureArea& rSrc); // no implementation given
	CSecureArea<T>& operator=(const CSecureArea& rSrc);  // no implementation given
};

template <class X> IRA::CSecAreaResourceWrapper<X>::CSecAreaResourceWrapper(BACIMutex *Mutex,X *Res): m_pLockMutex(Mutex), m_pLockResource(Res), m_bReleased(true) {
	m_pLockMutex->acquire();
	m_bReleased=false;
}

//copy constructor
template <class X> IRA::CSecAreaResourceWrapper<X>::CSecAreaResourceWrapper(const CSecAreaResourceWrapper& rSrc ) {
	m_pLockMutex=rSrc.m_pLockMutex;
	m_pLockResource=rSrc.m_pLockResource;
	m_bReleased=true;
	m_pLockMutex->acquire();
	m_bReleased=false;
}

template <class X> CSecAreaResourceWrapper<X>::~CSecAreaResourceWrapper() {
	Release();
}

template <class X> void CSecAreaResourceWrapper<X>::Release() {
	if (!m_bReleased) {
		m_pLockMutex->release();
		m_pLockResource=NULL;
		m_bReleased=true;
	}
}

template <class X> X* CSecAreaResourceWrapper<X>::operator->() {
	assert(m_pLockResource!=NULL);
	return m_pLockResource;
}

template <class X> X& CSecAreaResourceWrapper<X>::operator*() {
	assert(m_pLockResource!=NULL);
	return *m_pLockResource;
}

template <class T> CSecureArea<T>::CSecureArea(T* Obj): m_pResource(Obj), m_resourceOwner("")
{
}

template <class T> CSecureArea<T>::CSecureArea(bool Alloc): m_pResource(NULL), m_resourceOwner("")
{
	if (Alloc) {
		try {
			m_pResource=(T *)new T();
		}
		catch (std::bad_alloc& ex) {
			m_pResource=NULL;
		}
	}
}

template <class T> CSecureArea<T>::~CSecureArea()
{
	if (m_pResource) {
		m_Mutex.acquire(); //before closing the resource it is safe to acquire the mutex
		delete m_pResource;
		m_pResource=NULL;
	}
}

template <class T> void CSecureArea<T>::Init(T* Obj)
{
	if (!m_pResource) {
		m_pResource=Obj;
	}
}

template <class T> CSecAreaResourceWrapper<T> CSecureArea<T>::Get()
{
	return CSecAreaResourceWrapper<T>(&m_Mutex,m_pResource);
}

template <class T> CSecAreaResourceWrapper<T> CSecureArea<T>::Get(IRA::CString owner,bool show)
{
	if (show) {
		if (m_resourceOwner=="") {
			printf("%s is locking mutex, previous owner was none",(const char *)owner);
		}
		else {
			printf("%s is locking mutex, previous owner was %s",(const char *)owner,(const char *)m_resourceOwner);
		}
	}
	CSecAreaResourceWrapper<T> lock(&m_Mutex,m_pResource);
	m_resourceOwner=owner;
	if (show) printf ("......done!\n");
	return lock;
}


}

#endif /* _SECUREAREA_H */
