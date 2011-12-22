
#include <RWStepAP214_RWAutoDesignApprovalAssignment.ixx>
#include <StepAP214_HArray1OfAutoDesignGeneralOrgItem.hxx>
#include <StepAP214_AutoDesignGeneralOrgItem.hxx>
#include <StepBasic_Approval.hxx>


#include <Interface_EntityIterator.hxx>


#include <StepAP214_AutoDesignApprovalAssignment.hxx>


RWStepAP214_RWAutoDesignApprovalAssignment::RWStepAP214_RWAutoDesignApprovalAssignment () {}

void RWStepAP214_RWAutoDesignApprovalAssignment::ReadStep
	(const Handle(StepData_StepReaderData)& data,
	 const Standard_Integer num,
	 Handle(Interface_Check)& ach,
	 const Handle(StepAP214_AutoDesignApprovalAssignment)& ent) const
{


	// --- Number of Parameter Control ---

	if (!data->CheckNbParams(num,2,ach,"auto_design_approval_assignment")) return;

	// --- inherited field : assignedApproval ---

	Handle(StepBasic_Approval) aAssignedApproval;
	data->ReadEntity(num, 1,"assigned_approval", ach, STANDARD_TYPE(StepBasic_Approval), aAssignedApproval);

	// --- own field : items ---

	Handle(StepAP214_HArray1OfAutoDesignGeneralOrgItem) aItems;
	StepAP214_AutoDesignGeneralOrgItem aItemsItem;
	Standard_Integer nsub2;
	if (data->ReadSubList (num,2,"items",ach,nsub2)) {
	  Standard_Integer nb2 = data->NbParams(nsub2);
	  aItems = new StepAP214_HArray1OfAutoDesignGeneralOrgItem (1, nb2);
	  for (Standard_Integer i2 = 1; i2 <= nb2; i2 ++) {
	    Standard_Boolean stat2 = data->ReadEntity
	         (nsub2,i2,"items",ach,aItemsItem);
	    if (stat2) aItems->SetValue(i2,aItemsItem);
	  }
	}

	//--- Initialisation of the read entity ---


	ent->Init(aAssignedApproval, aItems);
}


void RWStepAP214_RWAutoDesignApprovalAssignment::WriteStep
	(StepData_StepWriter& SW,
	 const Handle(StepAP214_AutoDesignApprovalAssignment)& ent) const
{

	// --- inherited field assignedApproval ---

	SW.Send(ent->AssignedApproval());

	// --- own field : items ---

	SW.OpenSub();
	for (Standard_Integer i2 = 1;  i2 <= ent->NbItems();  i2 ++) {
	  SW.Send(ent->ItemsValue(i2).Value());
	}
	SW.CloseSub();
}


void RWStepAP214_RWAutoDesignApprovalAssignment::Share(const Handle(StepAP214_AutoDesignApprovalAssignment)& ent, Interface_EntityIterator& iter) const
{

	iter.GetOneItem(ent->AssignedApproval());


	Standard_Integer nbElem2 = ent->NbItems();
	for (Standard_Integer is2=1; is2<=nbElem2; is2 ++) {
	  iter.GetOneItem(ent->ItemsValue(is2).Value());
	}

}

