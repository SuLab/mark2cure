IntroStep = Backbone.Marionette.View.extend({
  // The initial page, classify + describe the activity
  template: '#intro-step',

  ui: {
    'next': '#next',
  },

  events: {
    'click @ui.next': function() { this.triggerMethod('step:complete:1', this); },
  },

});

Step1 = Backbone.Marionette.View.extend({
  // The initial page, classify + describe the activity
  template: '#step-one',

  ui: {
    'next': '#next',
  },

  regions: {
    'interactive': '#interactive'
  },

  events: {
    'click @ui.next': function() { this.triggerMethod('step:complete:2', this); },
  },

  onRender: function() {
    var main = this.getRegion('interactive');

    var collection = new REExtractionList([{
        "id": 0,
        "document": 0,
        "relation_type": "g_d",
        "concepts": {
          "c2": {
            "text": "Physics",
            "type": "d",
            "id": "0"
          },
          "c1": {
            "text": "Volunteers",
            "type": "g",
            "id": "0"
          }
        },
        "user_completed": false}]);
    var tree = new Tree({'collection': collection});
    main.show(tree.render());

    // var sub_view = new Tree({'collection': collection});
    // sub_view.render();

  //
  //   Tree['start'].onShow = function(m, c) {
  //     c.$el.css({'opacity': .25});
  //     c.$el.append('<div id="tree-cover"></div>');
  //
  //     #<{(| (TODO) Integrate into Tree app |)}>#
  //     var current_relationship = collection.findWhere({'current': true});
  //     add_relation_classes(current_relationship);
  //
  //     $el = $('#c1');
  //     $el.popover('hide');
  //     $el.popover({
  //       container: 'body',
  //       html: true,
  //       animation: false,
  //       content: function() {
  //           return 'A concept is a category to which a term can be assigned';
  //       },
  //       placement: 'top'
  //     });
  //     $el.popover('show');
  //
  //     $el2 = $('#c2');
  //     $el2.popover('hide');
  //     $el2.popover({
  //       container: 'body',
  //       html: true,
  //       animation: false,
  //       content: function() {
  //           return 'Eg- ‘Physics’ is part of the ‘field of science’ concept';
  //       },
  //       placement: 'top'
  //     });
  //     $el2.popover('show');
  //   }

  }
});


RelationTrainingOne = Backbone.Marionette.View.extend({
  el: '#relation-training-one',
  template: false,

  regions: {
    'one': '#step-1',
    'two': '#step-2',
    'three': '#step-3',
    'four': '#step-4',
    'five': '#step-5'
  },

  childViewEvents: {
    'step:complete:1': function(childView) {
      this.emptyRegions();
      this.showChildView('two', new Step1());
    },
    'step:complete:2': function(childView) {
      this.emptyRegions();
      this.showChildView('three', new Step2());
    },
  },

  onRender: function() {
    // Initial load of the page
      this.showChildView('one', new IntroStep());
  }

});

if($('#relation-training-one').length) {
  var sub_view = new RelationTrainingOne();
  sub_view.render();
}
