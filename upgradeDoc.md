This file documents the changes involved when I upgraded MarionetteJS from v2 to v3

# YPet

First, some changes have been made in `ypet.js`.

  * The initializor in ypet.js has been removed. Instead, initialize `YPet` before calling `YPet.start()` in each javascript file.
  * After initializing `YPet`, call `initAnnotationTypes(YPet)` to initialize the (default) annotation types. You don't need to call this if you are to manually set `YPet.AnnotationTypes` later.
  
## Initialize YPet

Since `addInitializer` has been removed in v3, I have created a same function as the callback function. For example,

```
YPet.addInitializer(function() {
  // Some code here...
});

YPet.start();
```

will be changed to

```
YPet = new Backbone.Marionette.Application({
  regions: "#someregion"
});

initAnnotationTypes(YPet);

var load_initial_data;
(load_initial_data = function() {
  // The same code in the callback of addInitializor
})();

YPet.start();
```

## Using `.getRegion` instead of `[]`

* If YPet has only one region, only a little change will be needed. For example, in `static/js/pages/training/entity-recognition/basic.js`,

  ```
  YPet.addInitializer(function() {
    // ...
    YPet.addRegions({'basics': '#basics'});
    YPet['basics'].show...
  });
  ```
  
  is changed to
  
  ```
  YPet = new Backbone.Marionette.Application({
    region: "#basics"
  });
  
  var load_initial_data;
  (load_initial_data = function() {
    // ...
    YPet.getRegion().show...
  })();
  ```
 
* If YPet has multiple regions, then v3 asks any app to have a view to hold multiple regions, so to update the following code, for example:

  ```
  YPet.addInitializer(function(options) {    
      var regions = {};

      _.each(passages, function(passage, passage_idx) {
        var passage_id = _.find(passage.infon, function(o){return o['@key']=='id';})['#text'];
        var p_body = '<div id="'+ passage_id +'" class="paragraph-box m-t-1"><p class="paragraph"></p></div></div>';
        $('.paragraphs').append(p_body);
        regions[''+passage_idx] = '#'+passage_id;
      });
      YPet.addRegions(regions);

      // To get a specific region
      YPet[passage_idx].show(...);
  );
  ```
  
  We notice that `.paragraphs` is where multiple editors are created, so it will be changed to
  
  ```
  YPet = new Backbone.Marionette.Application({
    region: "#doc_{{document.pk}}" // where the application is placed
  });
  initAnnotationTypes(YPet);
  
  var load_initial_data;
  (load_initial_data = function() {
    var regions = {};
    _.each(
      // ... same as above
    );
    
    // Instead of calling `YPet.addRegions(...)`, we add the view containing these regions to YPet
    YPet.showView(new (Mn.View.extend({
      el: ".paragraphs", // the parent of multiple editors, have to be specified
      regions: regions
    }))());
    
    // To get the region, instead of using `[]`, we use:
    YPet.getView().getRegion(passage_idx).show(...);
  })();
  ```
## Using `View` instead of `ItemView`

Simply replacing `ItemView` in the code with `View`. For example (in `js/libs/tree.js`), 

```
RelationView = Backbone.Marionette.ItemView.extend({ 
  ...
```
can be changed to
```
RelationView = Backbone.Marionette.View.extend({ 
  ...
```
