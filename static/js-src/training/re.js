var channel = Backbone.Radio.channel('tree');

//
// Views
//




RETrainingStepView = TrainingStepView.extend({

  onRender: function() {
    this.showChildView('text', new TrainingStepTextView({'model': this.model}));
    // this.showChildView('action', new RETrainingAction({collection: this.collection}));
    this.showChildView('footer', new TrainingFooterView({'model': this.model}));
  }

});


RETrainingModuleView = TrainingModuleView.extend({

  onRender: function() {
    this.showChildView('progress', new TrainingStepProgress({'collection': this.model.get('steps')}));
    var step = this.model.get('steps').get_active();
    this.showChildView('step', new RETrainingStepView({'model': step}));
  }

});


RETrainingTaskView = TrainingTaskView.extend({

  initialize: function() {
    this.collection = new TrainingModuleCollection([], {'task_type': 're'});
    this.collection.fetch();
  },

  onRender: function() {
    if(this.collection.length) {
      this.showChildView('navigation', new TrainingModuleNavigation({'collection': this.collection}));

      var module = this.collection.get_active();
      this.showChildView('module', new RETrainingModuleView({'model': module}));

    }
  }

});

