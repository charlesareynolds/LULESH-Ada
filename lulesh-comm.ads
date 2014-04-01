package LULESH.Comm is

   type Domain_Member is private;
   type Domain_Member_Array is array
     (Natural range <>) of Domain_Member;

   --- // lulesh-comm

   --x void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
   --x               Index_t dx, Index_t dy, Index_t dz,
   --x               bool doRecv, bool planeOnly);
   procedure Recv
     (domain    : in out Domain_Record;
      msgType    : in Int_t;
      xferFields : in Index_Type;
      dx         : in Element_Index;
      dy         : in Element_Index;
      dz         : in Element_Index;
      doRecv     : in Boolean;
      planeOnly  : in Boolean);

   --x void CommSend(Domain& domain, Int_t msgType,
   --x               Index_t xferFields, Domain_member *fieldData,
   --x               Index_t dx, Index_t dy, Index_t dz,
   --x               bool doSend, bool planeOnly);
   procedure Send
     (domain    : in out Domain_Record;
      msgType    : in Int_t;
      xferFields : in Index_Type;
      fieldData  : in Domain_Member_Array;
      dx         : in Element_Index;
      dy         : in Element_Index;
      dz         : in Element_Index;
      doSend     : in Boolean;
      planeOnly  : in Boolean);

   --x void CommSBN(Domain& domain, Int_t xferFields, Domain_member *fieldData);
   procedure SBN
     (domain    : in out Domain_Record;
      xferFields : in Int_t;
      fieldData  : out Domain_member);

   --x void CommSyncPosVel(Domain& domain);
   procedure SyncPosVel
     (domain : in out Domain_Record);

   --x void CommMonoQ(Domain& domain);
   procedure MonoQ
     (domain : in out Domain_Record);


private

   type Domain_Member is new Integer;

end LULESH.Comm;
